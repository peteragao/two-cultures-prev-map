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
X_pop <- readRDS(paste0("data/Zambia/zmb_X_pop.rds"))

# transform covariates
X_pop$l1pntl <- log1p(X_pop$ntl)
X_pop$l1paccess <- log1p(X_pop$access)
X_pop$urban <- 1 * (X_pop$urban)

# compute area-averaged covariates and weights by population
X_pop <- X_pop |>
  group_by(admin1_name) |>
  mutate(adm1_pop_weight = svy_f1549_pop / sum(svy_f1549_pop)) |>
  ungroup() |>
  group_by(admin2_name) |>
  mutate(adm2_pop_weight = svy_f1549_pop / sum(svy_f1549_pop)) |>
  ungroup()
X_pop$l1pntl_adm1 <- sum(X_pop$l1pntl * X_pop$adm1_pop_weight, na.rm = T)
X_pop$l1pntl_adm2 <- sum(X_pop$l1pntl * X_pop$adm2_pop_weight, na.rm = T)
X_pop$l1paccess_adm1 <- sum(X_pop$l1paccess * X_pop$adm1_pop_weight, na.rm = T)
X_pop$l1paccess_adm2 <- sum(X_pop$l1paccess * X_pop$adm2_pop_weight, na.rm = T)
X_pop$malaria_adm1 <- sum(X_pop$malaria * X_pop$adm1_pop_weight, na.rm = T)
X_pop$malaria_adm2 <- sum(X_pop$malaria * X_pop$adm2_pop_weight, na.rm = T)
X_adm1_avg <- X_pop |>
  group_by(admin1_name) |>
  mutate(l1pntl_adm1 = sum(l1pntl * adm1_pop_weight, na.rm = T),
         l1paccess_adm1 = sum(l1paccess * adm1_pop_weight, na.rm = T),
         malaria_adm1 = sum(malaria * adm1_pop_weight, na.rm = T)) |>
  select(admin1_name, l1pntl_adm1, l1paccess_adm1, malaria_adm1) |>
  slice(which.min(l1pntl_adm1)) |>
  ungroup() |>
  st_set_geometry(NULL)
X_adm2_avg <- X_pop |>
  group_by(admin2_name) |>
  mutate(l1pntl_adm2 = sum(l1pntl * adm2_pop_weight, na.rm = T),
         l1paccess_adm2 = sum(l1paccess * adm2_pop_weight, na.rm = T),
         malaria_adm2 = sum(malaria * adm2_pop_weight, na.rm = T)) |>
  select(admin2_name, l1pntl_adm2, l1paccess_adm2, malaria_adm2) |>
  slice(which.min(l1pntl_adm2)) |>
  ungroup() |>
  st_set_geometry(NULL)

# DHS RESPONSE DATA ------------------------------------------------------------
if (T) {
  svy_dat <- readRDS("data/Zambia/clean_DHS_data.rds")
} else {
  # Load cluster (EA) locations
  ea_locs <- st_read(geo_file)
  st_crs(ea_locs) <- st_crs(4326)
  # Remove EAs with missing geo information
  ea_locs <- st_crop(ea_locs, poly_adm0)
  
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
    stratum = ir_dat$v023
  )
  svy_dat <- merge(svy_dat, join_dat, by = c("cluster", "hshold", "line"))
  
  # merge geographic info 
  ea_locs <- ea_locs |> 
    select(DHSCLUST, URBAN_RURA) |>
    rename(cluster = DHSCLUST, urban = URBAN_RURA)
  
  svy_dat <- ea_locs |>
    right_join(svy_dat, by = "cluster") |>
    filter(!st_is_empty(geometry))
  
  # assign points to admin 2
  svy_dat <- st_join(svy_dat, poly_adm2 |> select(NAME_1, NAME_2),
                     join = st_nearest_feature)
  
  svy_dat <- svy_dat |>
    rename(admin1_name = NAME_1, admin2_name = NAME_2) |>
    mutate(admin1 = match(admin1_name, poly_adm1$NAME_1),
           admin2 = match(admin2_name, poly_adm2$NAME_2),
           admin1_char = paste0("admin1_", admin1),
           admin2_char = paste0("admin2_", admin2))
  
  svy_dat$access <- 
    get_cov("data/Zambia/Covariates/2015_accessibility_to_cities_v1.0.tif", svy_dat)
  
  svy_dat$ntl <- 
    get_cov("data/Zambia/Covariates/zmb_viirs_100m_2016.tif", svy_dat)
  
  svy_dat$malaria <- 
    get_cov("data/Zambia/Covariates/202206_Global_Pf_Incidence_Rate_ZMB_2018.tiff", svy_dat)
  
  svy_dat$urban_prob <- 
    get_cov("data/Zambia/Covariates/UR/urban_prob_surface.tif", svy_dat)
  
  svy_dat$urban = 1 * (svy_dat$urban == "U")
  svy_dat$l1pntl <- log1p(svy_dat$ntl)
  svy_dat$l1paccess <- log1p(svy_dat$access)
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
sample_des <- svydesign(id = ~cluster + hshold,
                        strata = ~stratum, nest = T,
                        weights = ~wt, data = svy_dat)
#### AREA LEVEL ADMIN1 ####
adm1_alm_res <- smoothArea(hiv ~ l1pntl_adm1 + l1paccess_adm1 + malaria_adm1,
                           domain = ~admin1_name, 
                           X.domain = X_adm1_avg,
                           design = sample_des, 
                           adj.mat = admin1_mat, 
                           transform = "logit", 
                           return.samples = T)
saveRDS(adm1_alm_res, "results/estimates/adm1_alm_res.rds")
#### AREA LEVEL ADMIN2 ####
adm2_alm_res <- smoothArea(hiv ~ l1pntl_adm2 + l1paccess_adm2 + malaria_adm2,
                           domain = ~admin2_name, 
                           design = sample_des, 
                           X.domain = X_adm2_avg,
                           adj.mat = admin2_mat, 
                           transform = "logit",
                           return.samples = T)
saveRDS(adm2_alm_res, "results/estimates/adm2_alm_res.rds")

#### AREA LEVEL ADMIN1 NO COVARIATES ####
adm1_alm_no_cov_res <- smoothArea(hiv ~ 1,
                                  domain = ~admin1_name, 
                                  X.domain = X_adm1_avg,
                                  design = sample_des, 
                                  adj.mat = admin1_mat, 
                                  transform = "logit", 
                                  return.samples = T)
saveRDS(adm1_alm_no_cov_res, "results/estimates/adm1_alm_no_cov_res.rds")

#### AREA LEVEL ADMIN2 NO COVARIATES ####
adm2_alm_no_cov_res <- smoothArea(hiv ~ 1,
                                  domain = ~admin2_name, 
                                  design = sample_des, 
                                  X.domain = X_adm2_avg,
                                  adj.mat = admin2_mat, 
                                  transform = "logit",
                                  return.samples = T)
saveRDS(adm2_alm_no_cov_res, "results/estimates/adm2_alm_no_cov_res.rds")

# BINOMIAL UNIT LEVEL MODELS ---------------------------------------------------

#### BINOMIAL ADMIN1 ####
bin_adm1_res <- 
  smoothUnit(
    hiv ~ urban + l1pntl*urban + l1paccess*urban + malaria*urban,
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
#### BINOMIAL ADMIN2 ####
bin_adm2_res <- 
  smoothUnit(
    hiv ~ urban + l1pntl*urban + l1paccess*urban + malaria*urban,
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

#### BINOMIAL ADMIN1 NO COVARIATES ####
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
saveRDS(bin_no_cov_adm1_res, "results/estimates/bin_no_cov_adm1_res.rds")
bin_no_cov_adm1_res$bym2.model.est$method <- "Binomial BYM2: No cov."
#### BINOMIAL ADMIN2 NO COVARIATES ####
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
bin_no_cov_adm2_res$bym2.model.est$method <- "Binomial BYM2: No cov."
saveRDS(bin_no_cov_adm2_res, "results/estimates/bin_no_cov_adm2_res.rds")
# LONOBINOMIAL UNIT LEVEL MODELS -----------------------------------------------
#### LONOBINOMIAL ADMIN1 ####
lono_adm1_res <- 
  fitLonobinomialBYM2(
    hiv ~ urban +  l1pntl*urban + l1paccess*urban + malaria*urban,
    domain = ~admin1_name, 
    design = sample_des, 
    adj.mat = admin1_mat,
    X.pop = X_pop,
    X.pop.weights = X_pop$adm1_pop_weight,
    return.samples = T
  )
lono_adm1_res$lonobinomial.bym2.model.est$method <- "Lonobinomial BYM2"
saveRDS(lono_adm1_res, "results/estimates/lono_adm1_res.rds")
#### LONOBINOMIAL ADMIN2 ####
lono_adm2_res <- 
  fitLonobinomialBYM2(
    hiv ~ urban +  l1pntl*urban + l1paccess*urban + malaria*urban,
    domain = ~admin2_name, 
    design = sample_des, 
    adj.mat = admin2_mat,
    X.pop = X_pop,
    X.pop.weights = X_pop$adm2_pop_weight,
    return.samples = T
  )
lono_adm2_res$lonobinomial.bym2.model.est$method <- "Lonobinomial BYM2"
saveRDS(lono_adm2_res, "results/estimates/lono_adm2_res.rds")

#### LONOBINOMIAL ADMIN1 NO COVARIATES ####
lono_no_cov_adm1_res <- 
  fitLonobinomialBYM2(
    hiv ~ urban,
    domain = ~admin1_name, 
    design = sample_des, 
    adj.mat = admin1_mat,
    X.pop = X_pop,
    X.pop.weights = X_pop$adm1_pop_weight,
    return.samples = T
  )
lono_no_cov_adm1_res$lonobinomial.bym2.model.est$method <- "Lonobinomial BYM2: No cov."
saveRDS(lono_no_cov_adm1_res, "results/estimates/lono_no_cov_adm1_res.rds")

#### LONOBINOMIAL ADMIN2 NO COVARIATES ####
lono_no_cov_adm2_res <- 
  fitLonobinomialBYM2(
    hiv ~ urban,
    domain = ~admin2_name, 
    design = sample_des, 
    adj.mat = admin2_mat,
    X.pop = X_pop,
    X.pop.weights = X_pop$adm2_pop_weight,
    return.samples = T
  )
lono_no_cov_adm2_res$lonobinomial.bym2.model.est$method <- "Lonobinomial BYM2: No cov."
saveRDS(lono_no_cov_adm2_res, "results/estimates/lono_no_cov_adm2_res.rds")

# BETABINOMIAL UNIT LEVEL MODELS -----------------------------------------------
#### BETABINOMIAL ADMIN1 ####
bbin_adm1_res <- 
  fitBetabinomialBYM2(
    hiv ~ urban +  l1pntl*urban + l1paccess*urban + malaria*urban,
    domain = ~admin1_name, 
    design = sample_des, 
    adj.mat = admin1_mat,
    X.pop = X_pop,
    X.pop.weights = X_pop$adm1_pop_weight,
    return.samples = T
  )
bbin_adm1_res$betabinomial.bym2.model.est$method <- "Betabinomial BYM2"
saveRDS(bbin_adm1_res, "results/estimates/bbin_adm1_res.rds")
#### BETABINOMIAL ADMIN2 ####
bbin_adm2_res <- 
  fitBetabinomialBYM2(
    hiv ~ urban +  l1pntl*urban + l1paccess*urban + malaria*urban,
    domain = ~admin2_name, 
    design = sample_des, 
    adj.mat = admin2_mat,
    X.pop = X_pop,
    X.pop.weights = X_pop$adm2_pop_weight,
    return.samples = T
  )
bbin_adm2_res$betabinomial.bym2.model.est$method <- "Betabinomial BYM2"
saveRDS(bbin_adm2_res, "results/estimates/bbin_adm2_res.rds")
#### BETABINOMIAL ADMIN1 NO COVARIATES ####
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
bbin_no_cov_adm1_res$betabinomial.bym2.model.est$method <- "Betabinomial BYM2: No cov."
saveRDS(bbin_no_cov_adm1_res, "results/estimates/bbin_no_cov_adm1_res.rds")
#### BETABINOMIAL ADMIN2 NO COVARIATES ####
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
bbin_no_cov_adm2_res$betabinomial.bym2.model.est$method <- "Betabinomial BYM2: No cov."
saveRDS(bbin_no_cov_adm2_res, "results/estimates/bbin_no_cov_adm2_res.rds")

# GEOSTATISTICAL UNIT MODELS ---------------------------------------------
sample_locs <- st_sample(poly_adm0, 1000)
mesh = inla.mesh.2d(loc.domain = st_coordinates(sample_locs),
                    max.edge = 0.25, offset = -0.1)


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
bbin_LGM_fit <- fitContLGM(formula = hiv ~ urban + l1pntl*urban + l1paccess*urban + malaria*urban,
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
                              X.pop = X_pop,
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

#### BETABINOMIAL SPDE NO COVARIATES ####
bbin_LGM_no_cov_fit <- 
  fitContLGM(formula = hiv ~ urban,
             family = "betabinomial",
             cluster = ~cluster,
             cluster.effect = F,
             data = svy_dat,
             mesh = mesh,
             pc.prior.range = prior.range,
             pc.prior.sigma = prior.sigma,
             pc.prior.clust = prior.clust)
saveRDS(bbin_LGM_no_cov_fit, "results/estimates/bbin_LGM_no_cov_fit.rds")
bbin_LGM_no_cov_res <- 
  smoothContLGM(bbin_LGM_no_cov_fit,
                X.pop = X_pop,
                domain = ~admin1_name + admin2_name,
                mesh,
                n.sample = 1000,
                cluster.effect = F,
                X.pop.weights = X_pop$adm1_pop_weight,
                level = 0.9,
                return.samples = T)
bbin_LGM_no_cov_res$binomial.spde.lgm.est$method <- "Betabinomial GRF: No cov."
bbin_LGM_no_cov_res$binomial.spde.lgm.est.subdomain$method <- "Betabinomial GRF: No cov."
saveRDS(bbin_LGM_no_cov_res, "results/estimates/bbin_LGM_no_cov_res.rds")

#### BINOMIAL SPDE ####
bin_LGM_fit <- 
  fitContLGM(formula = hiv ~ urban + l1pntl*urban + l1paccess*urban + malaria*urban,
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
                X.pop = X_pop,
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

#### BINOMIAL SPDE NO COVARIATES ####
bin_LGM_no_cov_fit <- 
  fitContLGM(formula = hiv ~ urban, 
             family = "binomial",
             cluster = ~cluster,
             cluster.effect = F,
             data = svy_dat,
             mesh = mesh,
             pc.prior.range = prior.range,
             pc.prior.sigma = prior.sigma,
             pc.prior.clust = prior.clust)
saveRDS(bin_LGM_no_cov_fit, "results/estimates/bin_LGM_no_cov_fit.rds")
bin_LGM_no_cov_res <- 
  smoothContLGM(bin_LGM_no_cov_fit,
                X.pop = X_pop,
                domain = ~admin1_name + admin2_name,
                mesh,
                n.sample = 1000,
                cluster.effect = F,
                X.pop.weights = X_pop$adm1_pop_weight,
                level = 0.9,
                return.samples = T)
bin_LGM_no_cov_res$binomial.spde.lgm.est$method <-
  "Binomial GRF: No cov."
bin_LGM_no_cov_res$binomial.spde.lgm.est.subdomain$method <-
  "Binomial GRF: No cov."
saveRDS(bin_LGM_no_cov_res, "results/estimates/bin_LGM_no_cov_res.rds")

#### LONOBINOMIAL SPDE ####
lono_LGM_fit <- 
  fitContLGM(formula = hiv ~ urban + l1pntl*urban + l1paccess*urban + malaria*urban,
             family = "binomial",
             cluster = ~cluster,
             cluster.effect = T,
             data = svy_dat,
             mesh = mesh,
             pc.prior.range = prior.range,
             pc.prior.sigma = prior.sigma,
             pc.prior.clust = prior.clust)

saveRDS(lono_LGM_fit, "results/estimates/lono_LGM_fit.rds")
lono_LGM_res <- 
  smoothContLGM(lono_LGM_fit,
                X.pop = X_pop,
                domain = ~admin1_name + admin2_name,
                mesh,
                n.sample = 1000,
                cluster.effect = T,
                X.pop.weights = X_pop$adm1_pop_weight,
                level = 0.9,
                return.samples = T)
lono_LGM_res$binomial.spde.lgm.est$method <- "Lonobinomial GRF"
lono_LGM_res$binomial.spde.lgm.est.subdomain$method <- "Lonobinomial GRF"
saveRDS(lono_LGM_res, "results/estimates/lono_LGM_res.rds")

#### LONOBINOMIAL SPDE NO COVARIATES ####
lono_LGM_no_cov_fit <- 
  fitContLGM(formula = hiv ~ urban, 
             family = "binomial",
             cluster = ~cluster,
             cluster.effect = T,
             data = svy_dat,
             mesh = mesh,
             pc.prior.range = prior.range,
             pc.prior.sigma = prior.sigma,
             pc.prior.clust = prior.clust)
saveRDS(lono_LGM_no_cov_fit, "results/estimates/lono_LGM_no_cov_fit.rds")
lono_LGM_no_cov_res <- 
  smoothContLGM(lono_LGM_no_cov_fit,
                X.pop = X_pop,
                domain = ~admin1_name + admin2_name,
                mesh,
                n.sample = 1000,
                cluster.effect = T,
                X.pop.weights = X_pop$adm1_pop_weight,
                level = 0.9,
                return.samples = T)
lono_LGM_no_cov_res$binomial.spde.lgm.est$method <-
  "Lonobinomial GRF: No cov."
lono_LGM_no_cov_res$binomial.spde.lgm.est.subdomain$method <-
  "Lonobinomial GRF: No cov."
saveRDS(lono_LGM_no_cov_res, "results/estimates/lono_LGM_no_cov_res.rds")

# GRID AGGREGATION -------------------------------------------------------------
easpa <- readr::read_csv("data/Zambia/DHS/zmb_2018_easpa.csv") |>
  rename(admin1_name = GADMarea) |>
  select(admin1_name, eaUrb, eaRur) |>
  pivot_longer(cols = c(eaUrb, eaRur), values_to = "n_ea") |>
  mutate(urban = 1 * (name == "eaUrb"))
easize <- readr::read_csv("data/Zambia/DHS/zmb_2018_easpa.csv") |>
  rename(admin1_name = GADMarea) |>
  select(admin1_name, avgSizeUrb, avgSizeRur) |>
  pivot_longer(cols = c(avgSizeUrb, avgSizeRur), values_to = "size") |>
  mutate(urban = 1 * (name == "avgSizeUrb"))
X_sim_frame <- simulateFrame(X_pop, easpa, easize)

#### LONOBINOMIAL ADMIN1 FRAME ####
lono_simframe_adm1_res <- 
  fitLonobinomialBYM2(
    hiv ~ urban +  l1pntl*urban + l1paccess*urban + malaria*urban,
    domain = ~admin1_name, 
    design = sample_des, 
    adj.mat = admin1_mat,
    X.pop = X_sim_frame,
    X.pop.weights = X_sim_frame$adm1_pop_weight,
    return.samples = T
  )
lono_simframe_adm1_res$lonobinomial.bym2.model.est$method <- "Lonobinomial BYM2 Sim. Frame"
saveRDS(lono_simframe_adm1_res, "results/estimates/lono_simframe_adm1_res.rds")
#### LONOBINOMIAL ADMIN2 FRAME ####
lono_simframe_adm2_res <- 
  fitLonobinomialBYM2(
    hiv ~ urban +  l1pntl*urban + l1paccess*urban + malaria*urban,
    domain = ~admin2_name, 
    design = sample_des, 
    adj.mat = admin2_mat,
    X.pop = X_sim_frame,
    X.pop.weights = X_sim_frame$adm2_pop_weight,
    return.samples = T
  )
lono_simframe_adm2_res$lonobinomial.bym2.model.est$method <- "Lonobinomial BYM2 Sim. Frame"
saveRDS(lono_simframe_adm2_res, "results/estimates/lono_simframe_adm2_res.rds")
#### BETABINOMIAL ADMIN1 FRAME ####
bbin_simframe_adm1_res <- 
  fitBetabinomialBYM2(
    hiv ~ urban +  l1pntl*urban + l1paccess*urban + malaria*urban,
    domain = ~admin1_name, 
    design = sample_des, 
    adj.mat = admin1_mat,
    X.pop = X_sim_frame,
    X.pop.weights = X_sim_frame$adm1_pop_weight,
    return.samples = T
  )
bbin_simframe_adm1_res$betabinomial.bym2.model.est$method <- "Betabinomial BYM2 Sim. Frame"
saveRDS(bbin_simframe_adm1_res, "results/estimates/bbin_simframe_adm1_res.rds")
#### BETABINOMIAL ADMIN2 FRAME ####
bbin_simframe_adm2_res <- 
  fitBetabinomialBYM2(
    hiv ~ urban +  l1pntl*urban + l1paccess*urban + malaria*urban,
    domain = ~admin2_name, 
    design = sample_des, 
    adj.mat = admin2_mat,
    X.pop = X_sim_frame,
    X.pop.weights = X_sim_frame$adm2_pop_weight,
    return.samples = T
  )
bbin_simframe_adm2_res$betabinomial.bym2.model.est$method <- "Betabinomial BYM2 Sim. Frame"
saveRDS(bbin_simframe_adm2_res, "results/estimates/bbin_simframe_adm2_res.rds")

#### LONOBINOMIAL SPDE FRAME ####
lono_LGM_simframe_res <- smoothContLGM(lono_LGM_fit,
                                       X.pop = X_sim_frame,
                                       domain = ~admin1_name + admin2_name,
                                       mesh,
                                       n.sample = 1000,
                                       cluster.effect = T,
                                       X.pop.weights = X_sim_frame$adm1_pop_weight,
                                       level = 0.9,
                                       return.samples = T)
lono_LGM_simframe_res$binomial.spde.lgm.est$method <- "Lonobinomial GRF Sim. Frame"
lono_LGM_simframe_res$binomial.spde.lgm.est.subdomain$method <- "Lonobinomial GRF Sim. Frame"
saveRDS(lono_LGM_simframe_res, "results/estimates/lono_LGM_simframe_res.rds")

#### LONOBINOMIAL SPDE FRAME NO COVARIATES ####
lono_LGM_no_cov_simframe_res <- smoothContLGM(lono_LGM_no_cov_fit,
                                              X.pop = X_sim_frame,
                                              domain = ~admin1_name + admin2_name,
                                              mesh,
                                              n.sample = 1000,
                                              cluster.effect = T,
                                              X.pop.weights = X_sim_frame$adm1_pop_weight,
                                              level = 0.9,
                                              return.samples = T)
lono_LGM_no_cov_simframe_res$binomial.spde.lgm.est$method <- "Lonobinomial GRF Sim. Frame: No cov."
lono_LGM_no_cov_simframe_res$binomial.spde.lgm.est.subdomain$method <- "Lonobinomial GRF Sim. Frame: No cov."
saveRDS(lono_LGM_no_cov_simframe_res, "results/estimates/lono_LGM_no_cov_simframe_res.rds")


#### BETABINOMIAL SPDE FRAME ####
bbin_LGM_simframe_res <- 
  smoothContLGM(bbin_LGM_fit,
                X.pop = X_sim_frame,
                domain = ~admin1_name + admin2_name,
                mesh,
                n.sample = 1000,
                cluster.effect = F,
                X.pop.weights = X_sim_frame$adm1_pop_weight,
                level = 0.9,
                return.samples = T)
bbin_LGM_simframe_res$betabinomial.spde.lgm.est$method <- 
  "Betabinomial GRF Sim. Frame"
bbin_LGM_simframe_res$betabinomial.spde.lgm.est.subdomain$method <- 
  "Betabinomial GRF Sim. Frame"
saveRDS(bbin_LGM_simframe_res, "results/estimates/bbin_LGM_simframe_res.rds")

#### BETABINOMIAL SPDE FRAME NO COVARIATES ####
bbin_LGM_no_cov_simframe_res <-
  smoothContLGM(bbin_LGM_no_cov_fit,
                X.pop = X_sim_frame,
                domain = ~admin1_name + admin2_name,
                mesh,
                n.sample = 1000,
                cluster.effect = F,
                X.pop.weights = X_sim_frame$adm1_pop_weight,
                level = 0.9,
                return.samples = T)
bbin_LGM_no_cov_simframe_res$betabinomial.spde.lgm.est$method <- 
  "Betabinomial GRF Sim. Frame: No cov."
bbin_LGM_no_cov_simframe_res$betabinomial.spde.lgm.est.subdomain$method <- 
  "Betabinomial GRF Sim. Frame: No cov."
saveRDS(bbin_LGM_no_cov_simframe_res, "results/estimates/bbin_LGM_no_cov_simframe_res.rds")


# FIGURES ----------------------------------------------------------------------
adm1_alm_res <- readRDS("results/estimates/adm1_alm_res.rds")
adm2_alm_res <- readRDS("results/estimates/adm2_alm_res.rds")
adm1_alm_no_cov_res <- readRDS("results/estimates/adm1_alm_no_cov_res.rds")
adm2_alm_no_cov_res <- readRDS("results/estimates/adm2_alm_no_cov_res.rds")

bin_adm1_res <- readRDS("results/estimates/bin_adm1_res.rds")
bin_adm2_res <- readRDS("results/estimates/bin_adm2_res.rds")
bin_no_cov_adm1_res <- readRDS("results/estimates/bin_no_cov_adm1_res.rds")
bin_no_cov_adm2_res <- readRDS("results/estimates/bin_no_cov_adm2_res.rds")

bbin_adm1_res <- readRDS("results/estimates/bbin_adm1_res.rds")
bbin_adm2_res <- readRDS("results/estimates/bbin_adm2_res.rds")
bbin_no_cov_adm1_res <- readRDS("results/estimates/bbin_no_cov_adm1_res.rds")
bbin_no_cov_adm2_res <- readRDS("results/estimates/bbin_no_cov_adm2_res.rds")
bbin_simframe_adm1_res <- readRDS("results/estimates/bbin_simframe_adm1_res.rds")
bbin_simframe_adm2_res <- readRDS("results/estimates/bbin_simframe_adm2_res.rds")

lono_adm1_res <- readRDS("results/estimates/lono_adm1_res.rds")
lono_adm2_res <- readRDS("results/estimates/lono_adm2_res.rds")
lono_no_cov_adm1_res <- readRDS("results/estimates/lono_no_cov_adm1_res.rds")
lono_no_cov_adm2_res <- readRDS("results/estimates/lono_no_cov_adm2_res.rds")
lono_simframe_adm1_res <- readRDS("results/estimates/lono_simframe_adm1_res.rds")
lono_simframe_adm2_res <- readRDS("results/estimates/lono_simframe_adm2_res.rds")

bin_LGM_res <- readRDS("results/estimates/bin_LGM_res.rds")
lono_LGM_res <- readRDS("results/estimates/lono_LGM_res.rds")
bbin_LGM_res <- readRDS("results/estimates/bbin_LGM_res.rds")
bin_LGM_no_cov_res <- readRDS("results/estimates/bin_LGM_no_cov_res.rds")
lono_LGM_no_cov_res <- readRDS("results/estimates/lono_LGM_no_cov_res.rds")
bbin_LGM_no_cov_res <- readRDS("results/estimates/bbin_LGM_no_cov_res.rds")
lono_LGM_simframe_res <- readRDS("results/estimates/lono_LGM_simframe_res.rds")
bbin_LGM_simframe_res <- readRDS("results/estimates/bbin_LGM_simframe_res.rds")

adm1_est_table <- bind_rows(
  adm1_alm_res$direct.est,
  adm1_alm_res$bym2.model.est |> mutate(method = "Fay-Herriot BYM2"),
  adm1_alm_no_cov_res$bym2.model.est |> mutate(method = "Fay-Herriot BYM2: No cov."),
  bin_adm1_res$bym2.model.est |> mutate(method = "Binomial BYM2"),
  bin_no_cov_adm1_res$bym2.model.est |> mutate(method = "Binomial BYM2: No cov."),
  lono_adm1_res$lonobinomial.bym2.model.est |> mutate(method = "Lonobinomial BYM2"),
  lono_no_cov_adm1_res$lonobinomial.bym2.model.est |> mutate(method = "Lonobinomial BYM2: No cov."),
  bbin_adm1_res$betabinomial.bym2.model.est |> mutate(method = "Betabinomial BYM2"),
  bbin_no_cov_adm1_res$betabinomial.bym2.model.est |> mutate(method = "Betabinomial BYM2: No cov."),
  bin_LGM_res$binomial.spde.lgm.est |> mutate(method = "Binomial GRF"),
  lono_LGM_res$binomial.spde.lgm.est |> mutate(method = "Lonobinomial GRF"),
  bbin_LGM_res$betabinomial.spde.lgm.est |> mutate(method = "Betabinomial GRF"),
  bin_LGM_no_cov_res$binomial.spde.lgm.est |> mutate(method = "Binomial GRF: No cov."),
  lono_LGM_no_cov_res$binomial.spde.lgm.est |> mutate(method = "Lonobinomial GRF: No cov."),
  bbin_LGM_no_cov_res$betabinomial.spde.lgm.est |> mutate(method = "Betabinomial GRF: No cov."),
  lono_simframe_adm1_res$lonobinomial.bym2.model.est |> mutate(method = "Lonobinomial BYM2 Sim Frame"),
  bbin_simframe_adm1_res$betabinomial.bym2.model.est |> mutate(method = "Betabinomial BYM2 Sim Frame"),
  lono_LGM_simframe_res$binomial.spde.lgm.est |> mutate(method = "Lonobinomial GRF Sim Frame"),
  bbin_LGM_simframe_res$betabinomial.spde.lgm.est |> mutate(method = "Lonobinomial GRF Sim Frame"),
) |>
  left_join(adm1_alm_res$direct.est |>
              select(domain, median, var) |>
              rename(direct.est = median, 
                     direct.est.var = var),
            by = "domain") |>
  mutate(CI = paste0("(", round(lower, 3), ", ", round(upper, 3), ")"))

write.csv(adm1_est_table, "results/estimates/adm1_est.csv")
write.csv(adm1_est_table, "internal/adm1_est.csv")

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

sample_des <- svydesign(id = ~cluster + hshold,
                        strata = ~stratum, nest = T,
                        weights = ~wt, data = svy_dat)

#### AREA LEVEL ADMIN1 ####
adm1_alm_res <- smoothArea(hiv ~ l1pntl_adm1 + l1paccess_adm1 + malaria_adm1,
                           domain = ~admin1_name, 
                           X.domain = X_adm1_avg,
                           design = sample_des, 
                           adj.mat = admin1_mat, 
                           transform = "logit", 
                           return.samples = T)
saveRDS(adm1_alm_res, "results/estimates/adm1_alm_res.rds")
#### AREA LEVEL ADMIN2 ####
adm2_alm_res <- smoothArea(hiv ~ l1pntl_adm2 + l1paccess_adm2 + malaria_adm2,
                           domain = ~admin2_name, 
                           design = sample_des, 
                           X.domain = X_adm2_avg,
                           adj.mat = admin2_mat, 
                           transform = "logit",
                           return.samples = T)
saveRDS(adm2_alm_res, "results/estimates/adm2_alm_res.rds")

#### AREA LEVEL ADMIN1 NO COVARIATES ####
adm1_alm_no_cov_res <- smoothArea(hiv ~ 1,
                                  domain = ~admin1_name, 
                                  X.domain = X_adm1_avg,
                                  design = sample_des, 
                                  adj.mat = admin1_mat, 
                                  transform = "logit", 
                                  return.samples = T)
saveRDS(adm1_alm_no_cov_res, "results/estimates/adm1_alm_no_cov_res.rds")

#### AREA LEVEL ADMIN2 NO COVARIATES ####
adm2_alm_no_cov_res <- smoothArea(hiv ~ 1,
                                  domain = ~admin2_name, 
                                  design = sample_des, 
                                  X.domain = X_adm2_avg,
                                  adj.mat = admin2_mat, 
                                  transform = "logit",
                                  return.samples = T)
saveRDS(adm2_alm_no_cov_res, "results/estimates/adm2_alm_no_cov_res.rds")

# BINOMIAL UNIT LEVEL MODELS ---------------------------------------------------

#### BINOMIAL ADMIN1 ####
bin_adm1_res <- 
  smoothUnit(
    hiv ~ urban + l1pntl*urban + l1paccess*urban + malaria*urban,
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
#### BINOMIAL ADMIN2 ####
bin_adm2_res <- 
  smoothUnit(
    hiv ~ urban + l1pntl*urban + l1paccess*urban + malaria*urban,
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

#### BINOMIAL ADMIN1 NO COVARIATES ####
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
saveRDS(bin_no_cov_adm1_res, "results/estimates/bin_no_cov_adm1_res.rds")
bin_no_cov_adm1_res$bym2.model.est$method <- "Binomial BYM2: No cov."
#### BINOMIAL ADMIN2 NO COVARIATES ####
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
bin_no_cov_adm2_res$bym2.model.est$method <- "Binomial BYM2: No cov."
saveRDS(bin_no_cov_adm2_res, "results/estimates/bin_no_cov_adm2_res.rds")
# LONOBINOMIAL UNIT LEVEL MODELS -----------------------------------------------
#### LONOBINOMIAL ADMIN1 ####
lono_adm1_res <- 
  fitLonobinomialBYM2(
    hiv ~ urban +  l1pntl*urban + l1paccess*urban + malaria*urban,
    domain = ~admin1_name, 
    design = sample_des, 
    adj.mat = admin1_mat,
    X.pop = X_pop,
    X.pop.weights = X_pop$adm1_pop_weight,
    return.samples = T
  )
lono_adm1_res$lonobinomial.bym2.model.est$method <- "Lonobinomial BYM2"
saveRDS(lono_adm1_res, "results/estimates/lono_adm1_res.rds")
#### LONOBINOMIAL ADMIN2 ####
lono_adm2_res <- 
  fitLonobinomialBYM2(
    hiv ~ urban +  l1pntl*urban + l1paccess*urban + malaria*urban,
    domain = ~admin2_name, 
    design = sample_des, 
    adj.mat = admin2_mat,
    X.pop = X_pop,
    X.pop.weights = X_pop$adm2_pop_weight,
    return.samples = T
  )
lono_adm2_res$lonobinomial.bym2.model.est$method <- "Lonobinomial BYM2"
saveRDS(lono_adm2_res, "results/estimates/lono_adm2_res.rds")

#### LONOBINOMIAL ADMIN1 NO COVARIATES ####
lono_no_cov_adm1_res <- 
  fitLonobinomialBYM2(
    hiv ~ urban,
    domain = ~admin1_name, 
    design = sample_des, 
    adj.mat = admin1_mat,
    X.pop = X_pop,
    X.pop.weights = X_pop$adm1_pop_weight,
    return.samples = T
  )
lono_no_cov_adm1_res$lonobinomial.bym2.model.est$method <- "Lonobinomial BYM2: No cov."
saveRDS(lono_no_cov_adm1_res, "results/estimates/lono_no_cov_adm1_res.rds")

#### LONOBINOMIAL ADMIN2 NO COVARIATES ####
lono_no_cov_adm2_res <- 
  fitLonobinomialBYM2(
    hiv ~ urban,
    domain = ~admin2_name, 
    design = sample_des, 
    adj.mat = admin2_mat,
    X.pop = X_pop,
    X.pop.weights = X_pop$adm2_pop_weight,
    return.samples = T
  )
lono_no_cov_adm2_res$lonobinomial.bym2.model.est$method <- "Lonobinomial BYM2: No cov."
saveRDS(lono_no_cov_adm2_res, "results/estimates/lono_no_cov_adm2_res.rds")

# BETABINOMIAL UNIT LEVEL MODELS -----------------------------------------------
#### BETABINOMIAL ADMIN1 ####
bbin_adm1_res <- 
  fitBetabinomialBYM2(
    hiv ~ urban +  l1pntl*urban + l1paccess*urban + malaria*urban,
    domain = ~admin1_name, 
    design = sample_des, 
    adj.mat = admin1_mat,
    X.pop = X_pop,
    X.pop.weights = X_pop$adm1_pop_weight,
    return.samples = T
  )
bbin_adm1_res$betabinomial.bym2.model.est$method <- "Betabinomial BYM2"
saveRDS(bbin_adm1_res, "results/estimates/bbin_adm1_res.rds")
#### BETABINOMIAL ADMIN2 ####
bbin_adm2_res <- 
  fitBetabinomialBYM2(
    hiv ~ urban +  l1pntl*urban + l1paccess*urban + malaria*urban,
    domain = ~admin2_name, 
    design = sample_des, 
    adj.mat = admin2_mat,
    X.pop = X_pop,
    X.pop.weights = X_pop$adm2_pop_weight,
    return.samples = T
  )
bbin_adm2_res$betabinomial.bym2.model.est$method <- "Betabinomial BYM2"
saveRDS(bbin_adm2_res, "results/estimates/bbin_adm2_res.rds")
#### BETABINOMIAL ADMIN1 NO COVARIATES ####
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
bbin_no_cov_adm1_res$betabinomial.bym2.model.est$method <- "Betabinomial BYM2: No cov."
saveRDS(bbin_no_cov_adm1_res, "results/estimates/bbin_no_cov_adm1_res.rds")
#### BETABINOMIAL ADMIN2 NO COVARIATES ####
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
bbin_no_cov_adm2_res$betabinomial.bym2.model.est$method <- "Betabinomial BYM2: No cov."
saveRDS(bbin_no_cov_adm2_res, "results/estimates/bbin_no_cov_adm2_res.rds")

# GEOSTATISTICAL UNIT MODELS ---------------------------------------------
sample_locs <- st_sample(poly_adm0, 1000)
mesh = inla.mesh.2d(loc.domain = st_coordinates(sample_locs),
                    max.edge = 0.25, offset = -0.1)


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
bbin_LGM_fit <- fitContLGM(formula = hiv ~ urban + l1pntl*urban + l1paccess*urban + malaria*urban,
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
                              X.pop = X_pop,
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

#### BETABINOMIAL SPDE NO COVARIATES ####
bbin_LGM_no_cov_fit <- 
  fitContLGM(formula = hiv ~ urban,
             family = "betabinomial",
             cluster = ~cluster,
             cluster.effect = F,
             data = svy_dat,
             mesh = mesh,
             pc.prior.range = prior.range,
             pc.prior.sigma = prior.sigma,
             pc.prior.clust = prior.clust)
saveRDS(bbin_LGM_no_cov_fit, "results/estimates/bbin_LGM_no_cov_fit.rds")
bbin_LGM_no_cov_res <- 
  smoothContLGM(bbin_LGM_no_cov_fit,
                X.pop = X_pop,
                domain = ~admin1_name + admin2_name,
                mesh,
                n.sample = 1000,
                cluster.effect = F,
                X.pop.weights = X_pop$adm1_pop_weight,
                level = 0.9,
                return.samples = T)
bbin_LGM_no_cov_res$binomial.spde.lgm.est$method <- "Betabinomial GRF: No cov."
bbin_LGM_no_cov_res$binomial.spde.lgm.est.subdomain$method <- "Betabinomial GRF: No cov."
saveRDS(bbin_LGM_no_cov_res, "results/estimates/bbin_LGM_no_cov_res.rds")

#### BINOMIAL SPDE ####
bin_LGM_fit <- 
  fitContLGM(formula = hiv ~ urban + l1pntl*urban + l1paccess*urban + malaria*urban,
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
                X.pop = X_pop,
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

#### BINOMIAL SPDE NO COVARIATES ####
bin_LGM_no_cov_fit <- 
  fitContLGM(formula = hiv ~ urban, 
             family = "binomial",
             cluster = ~cluster,
             cluster.effect = F,
             data = svy_dat,
             mesh = mesh,
             pc.prior.range = prior.range,
             pc.prior.sigma = prior.sigma,
             pc.prior.clust = prior.clust)
saveRDS(bin_LGM_no_cov_fit, "results/estimates/bin_LGM_no_cov_fit.rds")
bin_LGM_no_cov_res <- 
  smoothContLGM(bin_LGM_no_cov_fit,
                X.pop = X_pop,
                domain = ~admin1_name + admin2_name,
                mesh,
                n.sample = 1000,
                cluster.effect = F,
                X.pop.weights = X_pop$adm1_pop_weight,
                level = 0.9,
                return.samples = T)
bin_LGM_no_cov_res$binomial.spde.lgm.est$method <-
  "Binomial GRF: No cov."
bin_LGM_no_cov_res$binomial.spde.lgm.est.subdomain$method <-
  "Binomial GRF: No cov."
saveRDS(bin_LGM_no_cov_res, "results/estimates/bin_LGM_no_cov_res.rds")

#### LONOBINOMIAL SPDE ####
lono_LGM_fit <- 
  fitContLGM(formula = hiv ~ urban + l1pntl*urban + l1paccess*urban + malaria*urban,
             family = "binomial",
             cluster = ~cluster,
             cluster.effect = T,
             data = svy_dat,
             mesh = mesh,
             pc.prior.range = prior.range,
             pc.prior.sigma = prior.sigma,
             pc.prior.clust = prior.clust)

saveRDS(lono_LGM_fit, "results/estimates/lono_LGM_fit.rds")
lono_LGM_res <- 
  smoothContLGM(lono_LGM_fit,
                X.pop = X_pop,
                domain = ~admin1_name + admin2_name,
                mesh,
                n.sample = 1000,
                cluster.effect = T,
                X.pop.weights = X_pop$adm1_pop_weight,
                level = 0.9,
                return.samples = T)
lono_LGM_res$binomial.spde.lgm.est$method <- "Lonobinomial GRF"
lono_LGM_res$binomial.spde.lgm.est.subdomain$method <- "Lonobinomial GRF"
saveRDS(lono_LGM_res, "results/estimates/lono_LGM_res.rds")

#### LONOBINOMIAL SPDE NO COVARIATES ####
lono_LGM_no_cov_fit <- 
  fitContLGM(formula = hiv ~ urban, 
             family = "binomial",
             cluster = ~cluster,
             cluster.effect = T,
             data = svy_dat,
             mesh = mesh,
             pc.prior.range = prior.range,
             pc.prior.sigma = prior.sigma,
             pc.prior.clust = prior.clust)
saveRDS(lono_LGM_no_cov_fit, "results/estimates/lono_LGM_no_cov_fit.rds")
lono_LGM_no_cov_res <- 
  smoothContLGM(lono_LGM_no_cov_fit,
                X.pop = X_pop,
                domain = ~admin1_name + admin2_name,
                mesh,
                n.sample = 1000,
                cluster.effect = T,
                X.pop.weights = X_pop$adm1_pop_weight,
                level = 0.9,
                return.samples = T)
lono_LGM_no_cov_res$binomial.spde.lgm.est$method <-
  "Lonobinomial GRF: No cov."
lono_LGM_no_cov_res$binomial.spde.lgm.est.subdomain$method <-
  "Lonobinomial GRF: No cov."
saveRDS(lono_LGM_no_cov_res, "results/estimates/lono_LGM_no_cov_res.rds")

# GRID AGGREGATION -------------------------------------------------------------
easpa <- readr::read_csv("data/Zambia/DHS/zmb_2018_easpa.csv") |>
  rename(admin1_name = GADMarea) |>
  select(admin1_name, eaUrb, eaRur) |>
  pivot_longer(cols = c(eaUrb, eaRur), values_to = "n_ea") |>
  mutate(urban = 1 * (name == "eaUrb"))
easize <- readr::read_csv("data/Zambia/DHS/zmb_2018_easpa.csv") |>
  rename(admin1_name = GADMarea) |>
  select(admin1_name, avgSizeUrb, avgSizeRur) |>
  pivot_longer(cols = c(avgSizeUrb, avgSizeRur), values_to = "size") |>
  mutate(urban = 1 * (name == "avgSizeUrb"))
X_sim_frame <- simulateFrame(X_pop, easpa, easize)

#### LONOBINOMIAL ADMIN1 FRAME ####
lono_simframe_adm1_res <- 
  fitLonobinomialBYM2(
    hiv ~ urban +  l1pntl*urban + l1paccess*urban + malaria*urban,
    domain = ~admin1_name, 
    design = sample_des, 
    adj.mat = admin1_mat,
    X.pop = X_sim_frame,
    X.pop.weights = X_sim_frame$adm1_pop_weight,
    return.samples = T
  )
lono_simframe_adm1_res$lonobinomial.bym2.model.est$method <- "Lonobinomial BYM2 Sim. Frame"
saveRDS(lono_simframe_adm1_res, "results/estimates/lono_simframe_adm1_res.rds")
#### LONOBINOMIAL ADMIN2 FRAME ####
lono_simframe_adm2_res <- 
  fitLonobinomialBYM2(
    hiv ~ urban +  l1pntl*urban + l1paccess*urban + malaria*urban,
    domain = ~admin2_name, 
    design = sample_des, 
    adj.mat = admin2_mat,
    X.pop = X_sim_frame,
    X.pop.weights = X_sim_frame$adm2_pop_weight,
    return.samples = T
  )
lono_simframe_adm2_res$lonobinomial.bym2.model.est$method <- "Lonobinomial BYM2 Sim. Frame"
saveRDS(lono_simframe_adm2_res, "results/estimates/lono_simframe_adm2_res.rds")
#### BETABINOMIAL ADMIN1 FRAME ####
bbin_simframe_adm1_res <- 
  fitBetabinomialBYM2(
    hiv ~ urban +  l1pntl*urban + l1paccess*urban + malaria*urban,
    domain = ~admin1_name, 
    design = sample_des, 
    adj.mat = admin1_mat,
    X.pop = X_sim_frame,
    X.pop.weights = X_sim_frame$adm1_pop_weight,
    return.samples = T
  )
bbin_simframe_adm1_res$betabinomial.bym2.model.est$method <- "Betabinomial BYM2 Sim. Frame"
saveRDS(bbin_simframe_adm1_res, "results/estimates/bbin_simframe_adm1_res.rds")
#### BETABINOMIAL ADMIN2 FRAME ####
bbin_simframe_adm2_res <- 
  fitBetabinomialBYM2(
    hiv ~ urban +  l1pntl*urban + l1paccess*urban + malaria*urban,
    domain = ~admin2_name, 
    design = sample_des, 
    adj.mat = admin2_mat,
    X.pop = X_sim_frame,
    X.pop.weights = X_sim_frame$adm2_pop_weight,
    return.samples = T
  )
bbin_simframe_adm2_res$betabinomial.bym2.model.est$method <- "Betabinomial BYM2 Sim. Frame"
saveRDS(bbin_simframe_adm2_res, "results/estimates/bbin_simframe_adm2_res.rds")

#### LONOBINOMIAL SPDE FRAME ####
lono_LGM_simframe_res <- smoothContLGM(lono_LGM_fit,
                                       X.pop = X_sim_frame,
                                       domain = ~admin1_name + admin2_name,
                                       mesh,
                                       n.sample = 1000,
                                       cluster.effect = T,
                                       X.pop.weights = X_sim_frame$adm1_pop_weight,
                                       level = 0.9,
                                       return.samples = T)
lono_LGM_simframe_res$binomial.spde.lgm.est$method <- "Lonobinomial GRF Sim. Frame"
lono_LGM_simframe_res$binomial.spde.lgm.est.subdomain$method <- "Lonobinomial GRF Sim. Frame"
saveRDS(lono_LGM_simframe_res, "results/estimates/lono_LGM_simframe_res.rds")

#### LONOBINOMIAL SPDE FRAME NO COVARIATES ####
lono_LGM_no_cov_simframe_res <- smoothContLGM(lono_LGM_no_cov_fit,
                                              X.pop = X_sim_frame,
                                              domain = ~admin1_name + admin2_name,
                                              mesh,
                                              n.sample = 1000,
                                              cluster.effect = T,
                                              X.pop.weights = X_sim_frame$adm1_pop_weight,
                                              level = 0.9,
                                              return.samples = T)
lono_LGM_no_cov_simframe_res$binomial.spde.lgm.est$method <- "Lonobinomial GRF Sim. Frame: No cov."
lono_LGM_no_cov_simframe_res$binomial.spde.lgm.est.subdomain$method <- "Lonobinomial GRF Sim. Frame: No cov."
saveRDS(lono_LGM_no_cov_simframe_res, "results/estimates/lono_LGM_no_cov_simframe_res.rds")


#### BETABINOMIAL SPDE FRAME ####
bbin_LGM_simframe_res <- 
  smoothContLGM(bbin_LGM_fit,
                X.pop = X_sim_frame,
                domain = ~admin1_name + admin2_name,
                mesh,
                n.sample = 1000,
                cluster.effect = F,
                X.pop.weights = X_sim_frame$adm1_pop_weight,
                level = 0.9,
                return.samples = T)
bbin_LGM_simframe_res$betabinomial.spde.lgm.est$method <- 
  "Betabinomial GRF Sim. Frame"
bbin_LGM_simframe_res$betabinomial.spde.lgm.est.subdomain$method <- 
  "Betabinomial GRF Sim. Frame"
saveRDS(bbin_LGM_simframe_res, "results/estimates/bbin_LGM_simframe_res.rds")

#### BETABINOMIAL SPDE FRAME NO COVARIATES ####
bbin_LGM_no_cov_simframe_res <-
  smoothContLGM(bbin_LGM_no_cov_fit,
                X.pop = X_sim_frame,
                domain = ~admin1_name + admin2_name,
                mesh,
                n.sample = 1000,
                cluster.effect = F,
                X.pop.weights = X_sim_frame$adm1_pop_weight,
                level = 0.9,
                return.samples = T)
bbin_LGM_no_cov_simframe_res$betabinomial.spde.lgm.est$method <- 
  "Betabinomial GRF Sim. Frame: No cov."
bbin_LGM_no_cov_simframe_res$betabinomial.spde.lgm.est.subdomain$method <- 
  "Betabinomial GRF Sim. Frame: No cov."
saveRDS(bbin_LGM_no_cov_simframe_res, "results/estimates/bbin_LGM_no_cov_simframe_res.rds")


# NESTED MODEL -------------------------------------------------------------


admin_key <- as.data.frame(poly_adm2[, c("NAME_1", "NAME_2")])
admin_key <- admin_key[order(admin_key$NAME_2),]

# nested version of admin2_mat[]
admin2_mat_nested <- admin2_mat
for (i in 1:nrow(admin2_mat)) {
  admin2_mat_nested[i, which(poly_adm2$NAME_1 != poly_adm2$NAME_1[i])] <- 0
  if (sum(admin2_mat_nested[i, ]) > 0) {
    admin2_mat_nested[i,] <- admin2_mat_nested[i,] / sum(admin2_mat_nested[i,])
  }
}
X_pop$admin1_name[X_pop$admin1_name == "North-Western"] <- "NorthWestern"
sample_des$variables$admin1_name[sample_des$variables$admin1_name == "North-Western"] <- "NorthWestern"
bin_nested_adm2_res <- 
  smoothUnitNested(
    hiv ~ admin1_name + urban + l1pntl*urban + l1paccess*urban + malaria*urban,
    domain = ~admin2_name, 
    design = sample_des, 
    adj.mat = admin2_mat_nested,
    X.pop = X_pop,
    family = "binomial",
    X.pop.weights = X_pop$adm2_pop_weight,
    return.samples = T
  )
bin_nested_adm2_res$bym2.model.est$method <- "Binomial BYM2 Nested"
saveRDS(bin_nested_adm2_res, "results/estimates/bin_nested_adm2_res.rds")


bbin_nested_adm2_res <- 
  fitBetabinomialBYM2(
    hiv ~ admin1_name + urban + l1pntl*urban + l1paccess*urban + malaria*urban,
    domain = ~admin2_name, 
    design = sample_des, 
    adj.mat = admin2_mat_nested,
    X.pop = X_pop,
    X.pop.weights = X_pop$adm2_pop_weight,
    return.samples = T
  )
bbin_nested_adm2_res$betabinomial.bym2.model.est$method <- "Betabinomial BYM2 Nested"
saveRDS(bbin_nested_adm2_res, "results/estimates/bbin_nested_adm2_res.rds")

