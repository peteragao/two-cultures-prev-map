# get_final_estimates.R  
# LIBRARIES AND CONSTANTS ------------------------------------------------------
library(sf)
library(tidyr)
library(dplyr)
library(survey)
library(ggplot2)
library(spdep)
library(terra)
library(tidyterra)
library(readstata13)
library(SUMMER)
library(INLA)
library(sampling)
library(tidyverse)
geo_file <- "data/Zambia/DHS/ZMGE71FL/"
ir_file <- "data/Zambia/DHS/ZMIR71DT/ZMIR71FL.DTA"
hiv_file <- "data/Zambia/DHS/ZMAR71DT/ZMAR71FL.DTA"
gadm_abbrev <- "ZMB"
poly_path <- "data/Zambia/GADM/gadm41_ZMB_shp"

source("analysis/functions.R")

# GADM POLYGON DATA ------------------------------------------------------------

# GADM version 4.1 is used 
poly_layer_adm0 <- paste('gadm41', gadm_abbrev,
                         '0', sep = "_")
poly_layer_adm1 <- paste('gadm41', gadm_abbrev,
                         '1', sep = "_")
poly_layer_adm2 <- paste('gadm41', gadm_abbrev,
                         '2', sep = "_")

poly_adm0 <- st_read(dsn = poly_path,
                     layer = as.character(poly_layer_adm0)) 
poly_adm1 <-  st_read(dsn = poly_path,
                      layer = as.character(poly_layer_adm1)) 
poly_adm2 <-  st_read(dsn = poly_path,
                      layer = as.character(poly_layer_adm2)) 
st_crs(poly_adm0) <- st_crs(poly_adm1)  <- st_crs(poly_adm2)

admin1_mat <- poly2nb(poly_adm1)
admin1_mat <- nb2mat(admin1_mat, zero.policy = TRUE)

admin2_mat <- poly2nb(poly_adm2)
admin2_mat <- nb2mat(admin2_mat, zero.policy = TRUE)

poly_adm1$domain <- poly_adm1$NAME_1
poly_adm2$domain <- poly_adm2$NAME_2

rownames(admin1_mat) <- poly_adm1$NAME_1
rownames(admin2_mat) <- poly_adm2$NAME_2
# LOAD POPULATION SURFACE ------------------------------------------------------
load("data/Zambia/cov_info.RData")
load("data/Zambia/pop_info.RData")

X_pop <- zambia_cov$natl.grid |>
  rename(urban = strata, 
         admin1_name = admin1.name,
         admin2_name = admin2.name.full) |>
  mutate(admin2_name = str_split(admin2_name, '_', simplify = TRUE)[,2]) |>
  group_by(admin1_name) |>
  mutate(adm1_pop_weight = Population / sum(Population)) |>
  ungroup() |>
  group_by(admin2_name) |>
  mutate(adm2_pop_weight = Population / sum(Population)) |>
  ungroup()

# compute area-averaged covariates and weights by population
X_pop <- X_pop |>
  group_by(admin1_name) |>
  mutate(adm1_pop_weight = Population / sum(Population)) |>
  ungroup() |>
  group_by(admin2_name) |>
  mutate(adm2_pop_weight = Population / sum(Population)) |>
  ungroup()

X_adm1_avg <- X_pop |>
  group_by(admin1_name) |>
  mutate(l1pntl_adm1 = sum(l1pntl * adm1_pop_weight, na.rm = T),
         ntl_adm1 = sum(ntl * adm1_pop_weight, na.rm = T),
         l1paccess_adm1 = sum(l1paccess * adm1_pop_weight, na.rm = T),
         access_adm1 = sum(access * adm1_pop_weight, na.rm = T),
         malaria_adm1 = sum(malaria * adm1_pop_weight, na.rm = T)) |>
  dplyr::select(admin1_name, l1pntl_adm1, access_adm1, malaria_adm1) |>
  slice(which.min(l1pntl_adm1)) 
X_adm2_avg <- X_pop |>
  group_by(admin2_name) |>
  mutate(l1pntl_adm2 = sum(l1pntl * adm2_pop_weight, na.rm = T),
         ntl_adm2 = sum(l1pntl * adm2_pop_weight, na.rm = T),
         l1paccess_adm2 = sum(l1paccess * adm2_pop_weight, na.rm = T),
         access_adm2 = sum(access * adm2_pop_weight, na.rm = T),
         malaria_adm2 = sum(malaria * adm2_pop_weight, na.rm = T)) |>
  dplyr::select(admin2_name, l1pntl_adm2, access_adm2, malaria_adm2) |>
  slice(which.min(l1pntl_adm2)) 

X_pop <- X_pop |>
  left_join(X_adm1_avg, by = "admin1_name") |>
  left_join(X_adm2_avg, by = "admin2_name")

# DHS RESPONSE DATA ------------------------------------------------------------
if (T) {
  svy_dat <- readRDS("data/Zambia/clean_DHS_data.rds")
} else {
  # Load cluster (EA) locations
  ea_locs <- st_read(geo_file)
  st_crs(ea_locs) <- st_crs(4326)
  # Remove EAs with missing geo information
  ea_locs <- st_crop(ea_locs, poly_adm0)
  load("data/Zambia/cov_info.RData")
  # Load recode data
  ir_dat <- read.dta13(paste0(ir_file))
  hiv_dat <- read.dta13(paste0(hiv_file))
  
  svy_dat <- data.frame(
    cluster = hiv_dat$hivclust, 
    hshold = hiv_dat$hivnumb,
    line = hiv_dat$hivline,
    wt = hiv_dat$hiv05 / 1000000,
    hiv = 1 * (as.numeric(hiv_dat$hiv03) %in% 2:4)
  )
  join_dat <- data.frame(
    cluster = ir_dat$v001,
    hshold = ir_dat$v002,
    line = ir_dat$v003,
    stratum = ir_dat$v023,
    ir_wt = ir_dat$v005 / 1000000
  )
  svy_dat <- merge(svy_dat, join_dat, by = c("cluster", "hshold", "line"))
  
  # merge geographic info 
  ea_locs <- ea_locs |> 
    dplyr::select(DHSCLUST, URBAN_RURA) |>
    rename(cluster = DHSCLUST, urban = URBAN_RURA) 
  ea_locs$malaria <- get_cov('data/Zambia/Covariates/202206_Global_Pf_Incidence_Rate_ZMB_2018.tiff',
                             ea_locs, assign_na = T, standardize = T, l1p = F)
  ea_locs <- ea_locs |>
    left_join(cov.unit |> dplyr::select(-malaria), by = "cluster")
  
  svy_dat <- ea_locs |>
    right_join(svy_dat, by = "cluster") |>
    filter(!st_is_empty(geometry))

  svy_dat$urban <- 1 * (svy_dat$urban == "U")
  # assign points to admin 2
  svy_dat <- st_join(svy_dat, poly_adm2 |> dplyr::select(NAME_1, NAME_2),
                     join = st_nearest_feature)
  
  svy_dat <- svy_dat |>
    rename(admin1_name = NAME_1, admin2_name = NAME_2) |>
    mutate(admin1 = match(admin1_name, poly_adm1$NAME_1),
           admin2 = match(admin2_name, poly_adm2$NAME_2),
           admin1_char = paste0("admin1_", admin1),
           admin2_char = paste0("admin2_", admin2))
  
  svy_dat$n_trials <- 1
  
  svy_dat <- svy_dat |>
    left_join(X_adm1_avg, by = "admin1_name") |>
    left_join(X_adm2_avg, by = "admin2_name")
  
  saveRDS(svy_dat,
          file = paste0("data/Zambia/clean_DHS_data.rds"))
  jittered_locs <- st_jitter(ea_locs, factor = .01)
  jittered_locs <- jittered_locs |>
    filter(st_within(jittered_locs, st_as_sf(poly_adm0), sparse = F)) |>
    mutate(urban = ifelse(urban == "U", "Urban", "Rural"))
  
  # plot map of jittered locations
  zmb_map <- ggplot(data = st_as_sf(poly_adm2)) +
    geom_sf(lwd = .08, fill = NA) + 
    geom_sf(data = st_as_sf(poly_adm1), fill = NA, lwd = .66) + 
    geom_sf(data = jittered_locs,
            aes(color = urban),
            shape = 3, alpha = 1, size = .75) +
    scale_color_manual(values = c("mediumblue", "tomato"), name = NULL) + 
    guides(colour = guide_legend(override.aes = list(size = 4, stroke = 2))) +
    theme_bw() + guides(fill="none") +
    theme(legend.position="bottom",
          legend.text=element_text(size=16),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border=element_blank(),
          axis.line=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks=element_blank()) + 
    xlab("") + ylab("")
  ggsave(paste0("results/figures/Zambia_map.png"),
         width = 8, height = 7)
  zmb_map <- ggplot(data = st_as_sf(poly_adm2)) +
    geom_sf(lwd = .08, fill = NA) + 
    geom_sf(data = st_as_sf(poly_adm1), fill = NA, lwd = .66) + 
    geom_sf(data = jittered_locs,
            aes(color = urban),
            shape = 3, alpha = 1, size = .75) +
    geom_sf_label(data = st_as_sf(poly_adm1), aes(label = NAME_1)) +
    scale_color_manual(values = c("mediumblue", "tomato"), name = NULL) + 
    guides(colour = guide_legend(override.aes = list(size = 4, stroke = 2))) +
    theme_bw() + guides(fill="none") +
    theme(legend.position="bottom",
          legend.text=element_text(size=16),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border=element_blank(),
          axis.line=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks=element_blank()) + 
    xlab("") + ylab("")
  ggsave(paste0("results/figures/Zambia_map_labeled.png"),
         width = 8, height = 7)
}



# DIRECT ESTIMATION AND AREA LEVEL MODELS --------------------------------------

# model for unit-level covariate models
cov_formula <- hiv ~ urban + l1pntl + access + malaria
sample_des <- svydesign(id = ~cluster + hshold,
                        strata = ~stratum, nest = T,
                        weights = ~wt, data = svy_dat)

#### AREA LEVEL ADMIN1 ####
adm1_alm_res <- smoothArea(hiv ~ 1,
                           domain = ~admin1_name, 
                           X.domain = X_adm1_avg,
                           design = sample_des, 
                           adj.mat = admin1_mat, 
                           transform = "logit", 
                           return.samples = T)
adm1_cov_alm_res <- smoothArea(hiv ~ l1pntl_adm1 + access_adm1 + malaria_adm1,
                           domain = ~admin1_name, 
                           X.domain = X_adm1_avg,
                           design = sample_des, 
                           adj.mat = admin1_mat, 
                           transform = "logit", 
                           return.samples = T)
saveRDS(adm1_alm_res, "results/estimates/adm1_alm_res.rds")
saveRDS(adm1_cov_alm_res, "results/estimates/adm1_cov_alm_res.rds")
#### AREA LEVEL ADMIN2 ####
adm2_alm_res <- smoothArea(hiv ~ 1,
                           domain = ~admin2_name, 
                           design = sample_des, 
                           X.domain = X_adm2_avg,
                           adj.mat = admin2_mat, 
                           transform = "logit",
                           return.samples = T)
adm2_cov_alm_res <- smoothArea(hiv ~ l1pntl_adm2 + access_adm2 + malaria_adm2,
                           domain = ~admin2_name, 
                           design = sample_des, 
                           X.domain = X_adm2_avg,
                           adj.mat = admin2_mat, 
                           transform = "logit",
                           return.samples = T)
saveRDS(adm2_cov_alm_res, "results/estimates/adm2_cov_alm_res.rds")


# BINOMIAL UNIT LEVEL MODELS ---------------------------------------------------

#### BINOMIAL ADMIN1 ####
bin_adm1_res <- 
  smoothUnit(
    cov_formula, 
    domain = ~admin1_name, 
    design = sample_des, 
    adj.mat = admin1_mat,
    X.pop = X_pop,
    family = "binomial",
    X.pop.weights = X_pop$adm1_pop_weight,
    return.samples = T
  )
bin_adm1_res$bym2.model.est$method <- "Binomial BYM2"
saveRDS(bin_adm1_res, "results/estimates/bin_adm1_res.rds")
#### BINOMIAL ADMIN1 NO COV ####
bin_no_cov_adm1_res <- 
  smoothUnit(
    hiv ~ urban, 
    domain = ~admin1_name, 
    design = sample_des, 
    adj.mat = admin1_mat,
    X.pop = X_pop,
    family = "binomial",
    X.pop.weights = X_pop$adm1_pop_weight,
    return.samples = T
  )
bin_no_cov_adm1_res$bym2.model.est$method <- "Binomial BYM2 No Cov."
saveRDS(bin_no_cov_adm1_res, "results/estimates/bin_no_cov_adm1_res.rds")
#### BINOMIAL ADMIN2 ####
bin_adm2_res <- 
  smoothUnit(
    cov_formula,
    domain = ~admin2_name, 
    design = sample_des, 
    adj.mat = admin2_mat,
    X.pop = X_pop,
    family = "binomial",
    X.pop.weights = X_pop$adm2_pop_weight,
    return.samples = T
  )
bin_adm2_res$bym2.model.est$method <- "Binomial BYM2"
saveRDS(bin_adm2_res, "results/estimates/bin_adm2_res.rds")
#### BINOMIAL ADMIN2 NO COV ####
bin_no_cov_adm2_res <- 
  smoothUnit(
    hiv ~ urban, 
    domain = ~admin2_name, 
    design = sample_des, 
    adj.mat = admin2_mat,
    X.pop = X_pop,
    family = "binomial",
    X.pop.weights = X_pop$adm2_pop_weight,
    return.samples = T
  )
bin_no_cov_adm2_res$bym2.model.est$method <- "Binomial BYM2 No Cov."
saveRDS(bin_no_cov_adm2_res, "results/estimates/bin_no_cov_adm2_res.rds")

# BETABINOMIAL UNIT LEVEL MODELS -----------------------------------------------
#### BETABINOMIAL ADMIN1 ####
bbin_adm1_res <- 
  fitBetabinomialBYM2(
    cov_formula,
    domain = ~admin1_name, 
    design = sample_des, 
    adj.mat = admin1_mat,
    X.pop = X_pop,
    X.pop.weights = X_pop$adm1_pop_weight,
    return.samples = T
  )
bbin_adm1_res$betabinomial.bym2.model.est$method <- "Betabinomial BYM2"
saveRDS(bbin_adm1_res, "results/estimates/bbin_adm1_res.rds")
#### BETABINOMIAL ADMIN1 NO COV ####
bbin_no_cov_adm1_res <- 
  fitBetabinomialBYM2(
    hiv ~ urban,
    domain = ~admin1_name, 
    design = sample_des, 
    adj.mat = admin1_mat,
    X.pop = X_pop,
    X.pop.weights = X_pop$adm1_pop_weight,
    return.samples = T
  )
bbin_no_cov_adm1_res$betabinomial.bym2.model.est$method <- "Betabinomial BYM2 No Cov."
saveRDS(bbin_no_cov_adm1_res, "results/estimates/bbin_no_cov_adm1_res.rds")
#### BETABINOMIAL ADMIN2 ####
bbin_adm2_res <- 
  fitBetabinomialBYM2(
    cov_formula,
    domain = ~admin2_name, 
    design = sample_des, 
    adj.mat = admin2_mat,
    X.pop = X_pop,
    X.pop.weights = X_pop$adm2_pop_weight,
    return.samples = T
  )
bbin_adm2_res$betabinomial.bym2.model.est$method <- "Betabinomial BYM2"
saveRDS(bbin_adm2_res, "results/estimates/bbin_adm2_res.rds")
#### BETABINOMIAL ADMIN2 NO COV ####
bbin_no_cov_adm2_res <- 
  fitBetabinomialBYM2(
    hiv ~ urban,
    domain = ~admin2_name, 
    design = sample_des, 
    adj.mat = admin2_mat,
    X.pop = X_pop,
    X.pop.weights = X_pop$adm2_pop_weight,
    return.samples = T
  )
bbin_no_cov_adm2_res$betabinomial.bym2.model.est$method <- "Betabinomial BYM2 No Cov."
saveRDS(bbin_no_cov_adm2_res, "results/estimates/bbin_no_cov_adm2_res.rds")
# GEOSTATISTICAL UNIT MODELS ---------------------------------------------
sample_locs <- st_sample(poly_adm0, 1000)
mesh = inla.mesh.2d(loc.domain = st_coordinates(sample_locs),
                    max.edge = 0.25, offset = -0.1)
X_pop_sf <- st_as_sf(X_pop, coords = c("LONGNUM", "LATNUM"), crs = "EPSG:4326")

# Mesh
png('results/figures/Zambia_mesh.png', width = 1200, height = 800)
plot(mesh, asp = 1, main = "")
plot(poly_adm0, lwd = 3, add = TRUE)
plot(mesh, asp = 1, main = "", add = TRUE)
points(st_coordinates(jittered_locs) + # add jitter for privacy
         rnorm(n = length(st_coordinates(jittered_locs)), sd = .025), 
       col = "red", pch = 16, cex = .85)
dev.off()

# Priors
prior.range = c(3, 0.5)
prior.sigma = c(0.5, 0.5)
prior.clust  = list(prec = list(prior = "pc.prec",
                                param = c(1, 0.05)))

#### BETABINOMIAL SPDE ####
bbin_LGM_fit <- fitContLGM(formula = cov_formula,
                           family = "betabinomial",
                           cluster = ~cluster,
                           cluster.effect = F,
                           data = svy_dat,
                           mesh = mesh,
                           pc.prior.range = prior.range,
                           pc.prior.sigma = prior.sigma,
                           pc.prior.clust = prior.clust)
saveRDS(bbin_LGM_fit, "results/estimates/bbin_LGM_fit.rds")
bbin_LGM_res <- smoothContLGM(bbin_LGM_fit,
                              X.pop = X_pop_sf,
                              domain = ~admin1_name + admin2_name,
                              mesh,
                              n.sample = 1000,
                              cluster.effect = F,
                              X.pop.weights = X_pop$adm1_pop_weight,
                              level = 0.9,
                              return.samples = T)
bbin_LGM_res$betabinomial.spde.lgm.est$method <- "Betabinomial GRF"
bbin_LGM_res$betabinomial.spde.lgm.est.subdomain$method <- "Betabinomial GRF"
saveRDS(bbin_LGM_res, "results/estimates/bbin_LGM_res.rds")
#### BETABINOMIAL SPDE NO COV ####
bbin_no_cov_LGM_fit <- fitContLGM(formula = hiv ~ urban,
                                  family = "betabinomial",
                                  cluster = ~cluster,
                                  cluster.effect = F,
                                  data = svy_dat,
                                  mesh = mesh,
                                  pc.prior.range = prior.range,
                                  pc.prior.sigma = prior.sigma,
                                  pc.prior.clust = prior.clust)
saveRDS(bbin_no_cov_LGM_fit, "results/estimates/bbin_no_cov_LGM_fit.rds")
bbin_no_cov_LGM_res <- smoothContLGM(bbin_no_cov_LGM_fit,
                              X.pop = X_pop_sf,
                              domain = ~admin1_name + admin2_name,
                              mesh,
                              n.sample = 1000,
                              cluster.effect = F,
                              X.pop.weights = X_pop$adm1_pop_weight,
                              level = 0.9,
                              return.samples = T)
bbin_no_cov_LGM_res$betabinomial.spde.lgm.est$method <- "Betabinomial GRF No Cov."
bbin_no_cov_LGM_res$betabinomial.spde.lgm.est.subdomain$method <- "Betabinomial GRF No Cov."
saveRDS(bbin_no_cov_LGM_res, "results/estimates/bbin_no_cov_LGM_res.rds")
#### BINOMIAL SPDE ####
bin_LGM_fit <- 
  fitContLGM(formula = cov_formula,
             family = "binomial",
             cluster = ~cluster,
             cluster.effect = F,
             data = svy_dat,
             mesh = mesh,
             pc.prior.range = prior.range,
             pc.prior.sigma = prior.sigma,
             pc.prior.clust = prior.clust)

saveRDS(bin_LGM_fit, "results/estimates/bin_LGM_fit.rds")
bin_LGM_res <- 
  smoothContLGM(bin_LGM_fit,
                X.pop = X_pop_sf,
                domain = ~admin1_name + admin2_name,
                mesh,
                n.sample = 1000,
                cluster.effect = F,
                X.pop.weights = X_pop$adm1_pop_weight,
                level = 0.9,
                return.samples = T)
bin_LGM_res$binomial.spde.lgm.est$method <- "Binomial GRF"
bin_LGM_res$binomial.spde.lgm.est.subdomain$method <- "Binomial GRF"
saveRDS(bin_LGM_res, "results/estimates/bin_LGM_res.rds")
#### BINOMIAL SPDE NO COV ####
bin_no_cov_LGM_fit <- fitContLGM(formula = hiv ~ urban,
                                  family = "binomial",
                                  cluster = ~cluster,
                                  cluster.effect = F,
                                  data = svy_dat,
                                  mesh = mesh,
                                  pc.prior.range = prior.range,
                                  pc.prior.sigma = prior.sigma,
                                  pc.prior.clust = prior.clust)
saveRDS(bin_no_cov_LGM_fit, "results/estimates/bin_no_cov_LGM_fit.rds")
bin_no_cov_LGM_res <- smoothContLGM(bin_no_cov_LGM_fit,
                                     X.pop = X_pop_sf,
                                     domain = ~admin1_name + admin2_name,
                                     mesh,
                                     n.sample = 1000,
                                     cluster.effect = F,
                                     X.pop.weights = X_pop$adm1_pop_weight,
                                     level = 0.9,
                                     return.samples = T)
bin_no_cov_LGM_res$binomial.spde.lgm.est$method <- "Binomial GRF No Cov."
bin_no_cov_LGM_res$binomial.spde.lgm.est.subdomain$method <- "Binomial GRF No Cov."
saveRDS(bin_no_cov_LGM_res, "results/estimates/bin_no_cov_LGM_res.rds")
# SUMMARIZE RESULTS ------------------------------------------------------------

adm1_alm_res <- readRDS("results/estimates/adm1_alm_res.rds")
adm2_alm_res <- readRDS("results/estimates/adm2_alm_res.rds")

bin_adm1_res <- readRDS("results/estimates/bin_adm1_res.rds")
bin_adm2_res <- readRDS("results/estimates/bin_adm2_res.rds")

bbin_adm1_res <- readRDS("results/estimates/bbin_adm1_res.rds")
bbin_adm2_res <- readRDS("results/estimates/bbin_adm2_res.rds")
bbin_no_cov_adm1_res <- readRDS("results/estimates/bbin_no_cov_adm1_res.rds")
bbin_no_cov_adm2_res <- readRDS("results/estimates/bbin_no_cov_adm2_res.rds")

bin_LGM_res <- readRDS("results/estimates/bin_LGM_res.rds")
bbin_LGM_res <- readRDS("results/estimates/bbin_LGM_res.rds")

bbin_LGM_no_cov_res <- readRDS("results/estimates/bbin_LGM_no_cov_res.rds")


adm1_est_table <- bind_rows(
  adm1_alm_res$direct.est,
  adm1_alm_res$iid.model.est |> mutate(method = "Fay-Herriot IID"),
  adm1_alm_res$bym2.model.est |> mutate(method = "Fay-Herriot BYM2"),
  bin_adm1_res$bym2.model.est |> mutate(method = "Binomial BYM2"),
  bbin_adm1_res$betabinomial.bym2.model.est |> mutate(method = "Betabinomial BYM2"),
  bbin_no_cov_adm1_res$betabinomial.bym2.model.est |> mutate(method = "Betabinomial BYM2 No Cov."),
  bin_LGM_res$binomial.spde.lgm.est |> mutate(method = "Binomial GRF"),
  bbin_LGM_res$betabinomial.spde.lgm.est |> mutate(method = "Betabinomial GRF"),
  bbin_LGM_no_cov_res$betabinomial.spde.lgm.est |> mutate(method = "Betabinomial GRF No Cov."),
) |>
  left_join(adm1_alm_res$direct.est |>
              select(domain, median, var) |>
              rename(direct.est = median, 
                     direct.est.var = var),
            by = "domain") |>
  mutate(CI = paste0("(", round(lower, 3), ", ", round(upper, 3), ")"))
write.csv(adm1_est_table, "results/estimates/adm1_est.csv")

adm2_est_table <- bind_rows(
  adm2_alm_res$direct.est,
  adm2_alm_res$iid.model.est |> mutate(method = "Fay-Herriot IID"),
  adm2_alm_res$bym2.model.est |> mutate(method = "Fay-Herriot BYM2"),
  bin_adm2_res$bym2.model.est |> mutate(method = "Binomial BYM2"),
  bbin_adm2_res$betabinomial.bym2.model.est |> mutate(method = "Betabinomial BYM2"),
  
  bbin_no_cov_adm2_res$betabinomial.bym2.model.est |> mutate(method = "Betabinomial BYM2 No Cov."),
  bin_LGM_res$binomial.spde.lgm.est.subdomain |> mutate(method = "Binomial GRF"),
  bbin_LGM_res$betabinomial.spde.lgm.est.subdomain |> mutate(method = "Betabinomial GRF"),
  bbin_LGM_no_cov_res$betabinomial.spde.lgm.est.subdomain |> mutate(method = "Betabinomial GRF No Cov.")
) |>
  left_join(adm2_alm_res$direct.est |>
              select(domain, median, var) |>
              rename(direct.est = median, 
                     direct.est.var = var),
            by = "domain") |>
  mutate(CI = paste0("(", round(lower, 3), ", ", round(upper, 3), ")"))
write.csv(adm2_est_table, "results/estimates/adm2_est.csv")



