#-----------------------------------------------------------------------------#
# admin2_cross_validation.R
# This                                                                          
#
# Input files:
# Output files:
# Reference: 
#   - https://dhsprogram.com/data/Guide-to-DHS-Statistics/index.htm#t=HIV_Prevalence.htm%23Among_women_and_menbc-1&rhtocid=_17_1_0
#------------------------------------------------------------------------------#

library(sf)
library(tidyr)
library(dplyr)
library(survey)
library(ggplot2)
library(spdep)
library(terra)
library(readstata13)
library(SUMMER)
library(INLA)
library(sampling)
library(stringr)
geo_file <- "data/Zambia/DHS/ZMGE71FL/"
ir_file <- "data/Zambia/DHS/ZMIR71DT/ZMIR71FL.DTA"
hiv_file <- "data/Zambia/DHS/ZMAR71DT/ZMAR71FL.DTA"
#source("analysis/functions.R")
source("~/prevalence-mapping-review/analysis/functions.R")
setwd("~/prevalence-mapping-review/")
#INLA:::inla.binary.install()
inla.setOption("num.threads", 8)
# GADM POLYGON DATA ------------------------------------------------------------
gadm_abbrev <- "ZMB"
poly_path <- "data/Zambia/GADM/gadm41_ZMB_shp"
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

X_pop_sf <- st_as_sf(X_pop, coords = c("LONGNUM", "LATNUM"), crs = "EPSG:4326")
# DHS RESPONSE DATA ------------------------------------------------------------
svy_dat <- readRDS("data/Zambia/clean_DHS_data.rds")

# ID HOLDOUT AREAS -------------------------------------------------------------
# not all admin-2 areas have reliable direct estimates
# keep those with direct estimates with CV < 20% and variance > 0
adm2_est_table <- read.csv("results/estimates/adm2_est.csv")
  

holdout_direct_est <- adm2_est_table |>
  dplyr::filter(method == "Direct") |> 
  dplyr::filter(!is.na(mean)) |>
  mutate(CV = sqrt(var) / median * 100) |>
  filter(!(domain %in% c("Chama", "Chipata", "Lunga", "Mpika", "Mpulungu")) & var > 1e-10)
select_levels <-
  holdout_direct_est$domain[order(holdout_direct_est$mean, decreasing = T)] 


# CROSS-VALIDATION -------------------------------------------------------------
full_res_list <- list()
options(survey.lonely.psu="adjust")

sample_locs <- st_sample(poly_adm0, 1000)
mesh = inla.mesh.2d(loc.domain = st_coordinates(sample_locs),
                    max.edge = 0.25, offset = -0.1)
prior.range = c(3, 0.5)
prior.sigma = c(0.5, 0.5)
prior.clust  = list(prec = list(prior = "pc.prec",
                                param = c(1, 0.05)))


cov_formula <- hiv ~ urban + l1pntl + access + malaria
m <- nrow(holdout_direct_est)
for (i in 1:m) {
  res_list <- list()
  holdout_name <- holdout_direct_est$domain[i]
  
  if (!file.exists(paste0("results/estimates/holdout_res_", holdout_name, ".rds"))) {
    
    cat("\nTake out region: ")
    cat(holdout_name, sep = "")
    cat("\n")
    print(Sys.time())
    
    holdout_dat <- svy_dat |>
      filter(admin2_name != holdout_name)
    holdout_des <- svydesign(id = ~cluster + hshold,
                             strata = ~stratum, nest = T,
                             weights = ~wt, data = holdout_dat)
    #### MODELS ####
    
    #### AREA LEVEL ####
    adm2_alm_res <- 
      smoothArea(hiv ~ l1pntl_adm2 + access_adm2 + malaria_adm2,
                 domain = ~admin2_name, 
                 X.domain = X_adm2_avg,
                 design = holdout_des, 
                 adj.mat = admin2_mat, 
                 transform = "logit", 
                 level = .8,
                 return.samples = T)
    
    holdout_i <- match(holdout_name, adm2_alm_res$iid.model.est$domain)
    iid_logit_samples <-
      SUMMER::logit(adm2_alm_res$iid.model.sample)
    adm2_alm_res$iid.model.est$logit_mean <- 
      apply(iid_logit_samples, 1, mean)
    adm2_alm_res$iid.model.est$logit_var <- 
      apply(iid_logit_samples, 1, var)
    adm2_alm_res$iid.model.est$method <- "Fay-Herriot IID"
    bym2_logit_samples <-
      SUMMER::logit(adm2_alm_res$bym2.model.sample)[match(unique(adm2_est_table$domain), rownames(admin2_mat)), ]
    adm2_alm_res$bym2.model.est$logit_mean <- 
      apply(bym2_logit_samples, 1, mean)
    adm2_alm_res$bym2.model.est$logit_var <- 
      apply(bym2_logit_samples, 1, var)
    adm2_alm_res$bym2.model.est$method <- "Fay-Herriot BYM2"
    res_list <- 
      c(res_list,
        list(adm2_alm_res$iid.model.est[holdout_i, ]),
        list(adm2_alm_res$bym2.model.est[holdout_i, ]))
    
    #### AREA LEVEL NO COVARIATES ####
    adm2_alm_no_cov_res <- 
      smoothArea(hiv~1, domain = ~admin2_name, 
                 design = holdout_des, 
                 adj.mat = admin2_mat, 
                 transform = "logit", 
                 level = .8,
                 return.samples = T)
    iid_logit_samples <-
      SUMMER::logit(adm2_alm_no_cov_res$iid.model.sample)
    adm2_alm_no_cov_res$iid.model.est$logit_mean <- 
      apply(iid_logit_samples, 1, mean)
    adm2_alm_no_cov_res$iid.model.est$logit_var <- 
      apply(iid_logit_samples, 1, var)
    adm2_alm_no_cov_res$iid.model.est$method <- "Fay-Herriot IID: No cov."
    bym2_logit_samples <-
      SUMMER::logit(adm2_alm_no_cov_res$bym2.model.sample)[match(unique(adm2_est_table$domain), rownames(admin2_mat)), ]
    adm2_alm_no_cov_res$bym2.model.est$logit_mean <- 
      apply(bym2_logit_samples, 1, mean)
    adm2_alm_no_cov_res$bym2.model.est$logit_var <- 
      apply(bym2_logit_samples, 1, var)
    adm2_alm_no_cov_res$bym2.model.est$method <- "Fay-Herriot BYM2: No cov."
    res_list <- 
      c(res_list,
        list(adm2_alm_no_cov_res$iid.model.est[holdout_i, ]),
        list(adm2_alm_no_cov_res$bym2.model.est[holdout_i, ]))
    
    #### BINOMIAL ####
    bin_adm2_res <-
      smoothUnit(
        cov_formula,
        domain = ~admin2_name, 
        design = holdout_des, 
        adj.mat = admin2_mat,
        X.pop = X_pop,
        family = "binomial",
        X.pop.weights = X_pop$adm2_pop_weight,
        level = .8,
        return.samples = T
      )
    bin_adm2_res$bym2.model.est$method <- "Binomial BYM2"
    holdout_i <- match(holdout_name, bin_adm2_res$bym2.model.est$domain)
    bym2_logit_samples <-
      SUMMER::logit(bin_adm2_res$bym2.model.sample)[match(unique(adm2_est_table$domain), rownames(admin2_mat)), ]
    bin_adm2_res$bym2.model.est$logit_mean <- 
      apply(bym2_logit_samples, 1, mean)
    bin_adm2_res$bym2.model.est$logit_var <- 
      apply(bym2_logit_samples, 1, var)
    res_list <- 
      c(res_list,
        list(bin_adm2_res$bym2.model.est[holdout_i, ]))
    
    #### BINOMIAL NO COVARIATES ####
    bin_no_cov_adm2_res <-
      smoothUnit(
        hiv ~ urban,
        domain = ~admin2_name, 
        design = holdout_des, 
        adj.mat = admin2_mat,
        X.pop = X_pop,
        family = "binomial",
        X.pop.weights = X_pop$adm2_pop_weight,
        level = .8,
        return.samples = T
      )
    bin_no_cov_adm2_res$bym2.model.est$method <- "Binomial BYM2: No cov."
    holdout_i <- match(holdout_name, bin_no_cov_adm2_res$bym2.model.est$domain)
    bym2_logit_samples <-
      SUMMER::logit(bin_no_cov_adm2_res$bym2.model.sample)[match(unique(adm2_est_table$domain), rownames(admin2_mat)), ]
    bin_no_cov_adm2_res$bym2.model.est$logit_mean <- 
      apply(bym2_logit_samples, 1, mean)
    bin_no_cov_adm2_res$bym2.model.est$logit_var <- 
      apply(bym2_logit_samples, 1, var)
    res_list <- 
      c(res_list,
        list(bin_no_cov_adm2_res$bym2.model.est[holdout_i, ]))
    
    
    #### BETABINOMIAL ####
    bbin_adm2_res <-
      fitBetabinomialBYM2(
        cov_formula,
        domain = ~admin2_name, 
        design = holdout_des, 
        adj.mat = admin2_mat,
        X.pop = X_pop,
        X.pop.weights = X_pop$adm2_pop_weight,
        level = .8,
        return.samples = T
      )
    bbin_adm2_res$betabinomial.bym2.model.est$method <- "Betabinomial BYM2"
    holdout_i <- match(holdout_name, bbin_adm2_res$betabinomial.bym2.model.est$domain)
    bym2_logit_samples <-
      SUMMER::logit(bbin_adm2_res$betabinomial.bym2.model.sample)[match(unique(adm2_est_table$domain), rownames(admin2_mat)), ]
    bbin_adm2_res$betabinomial.bym2.model.est$logit_mean <- 
      apply(bym2_logit_samples, 1, mean)
    bbin_adm2_res$betabinomial.bym2.model.est$logit_var <- 
      apply(bym2_logit_samples, 1, var)
    res_list <- 
      c(res_list,
        list(bbin_adm2_res$betabinomial.bym2.model.est[holdout_i, ]))
    #### BETABINOMIAL NO COVARIATES ####
    bbin_no_cov_adm2_res <-
      fitBetabinomialBYM2(
        hiv ~ urban,
        domain = ~admin2_name, 
        design = holdout_des, 
        adj.mat = admin2_mat,
        X.pop = X_pop,
        X.pop.weights = X_pop$adm2_pop_weight,
        level = .8,
        return.samples = T
      )
    bbin_no_cov_adm2_res$betabinomial.bym2.model.est$method <-
      "Betabinomial BYM2: No cov."
    holdout_i <- match(holdout_name, bbin_no_cov_adm2_res$betabinomial.bym2.model.est$domain)
    bym2_logit_samples <-
      SUMMER::logit(bbin_no_cov_adm2_res$betabinomial.bym2.model.sample)[match(unique(adm2_est_table$domain), rownames(admin2_mat)), ]
    bbin_no_cov_adm2_res$betabinomial.bym2.model.est$logit_mean <- 
      apply(bym2_logit_samples, 1, mean)
    bbin_no_cov_adm2_res$betabinomial.bym2.model.est$logit_var <- 
      apply(bym2_logit_samples, 1, var)
    res_list <- 
      c(res_list,
        list(bbin_no_cov_adm2_res$betabinomial.bym2.model.est[holdout_i, ]))
    
    #### BINOMIAL SPDE ####
    bin_LGM_fit <- 
      fitContLGM(formula = cov_formula,
                 family = "binomial",
                 cluster = ~cluster,
                 cluster.effect = F,
                 data = holdout_dat,
                 mesh = mesh,
                 pc.prior.range = prior.range,
                 pc.prior.sigma = prior.sigma,
                 pc.prior.clust = prior.clust)
    bin_LGM_res <- 
      smoothContLGM(bin_LGM_fit,
                    X.pop = X_pop_sf,
                    domain = ~admin1_name + admin2_name,
                    mesh,
                    n.sample = 1000,
                    cluster.effect = F,
                    X.pop.weights = X_pop$adm2_pop_weight,
                    level = .8,
                    return.samples = T)
    bin_LGM_res$binomial.spde.lgm.est.subdomain$method <- "Binomial GRF"
    holdout_i <- match(holdout_name,  bin_LGM_res$binomial.spde.lgm.est.subdomain$domain)
    
    spde_logit_samples <-
      SUMMER::logit(bin_LGM_res$binomial.spde.lgm.sample.subdomain)
    bin_LGM_res$binomial.spde.lgm.est.subdomain$logit_mean <- 
      apply(spde_logit_samples, 1, mean)
    bin_LGM_res$binomial.spde.lgm.est.subdomain$logit_var <- 
      apply(spde_logit_samples, 1, var)

    res_list <- 
      c(res_list,
        list(bin_LGM_res$binomial.spde.lgm.est.subdomain[holdout_i, ]))

   
    #### BETABINOMIAL SPDE ####
    print("GRF MODELS")
    bbin_LGM_fit <- 
      fitContLGM(formula = cov_formula,
                 family = "betabinomial",
                 cluster = ~cluster,
                 cluster.effect = F,
                 data = holdout_dat,
                 mesh = mesh,
                 pc.prior.range = prior.range,
                 pc.prior.sigma = prior.sigma,
                 pc.prior.clust = prior.clust)
    bbin_LGM_res <- 
      smoothContLGM(bbin_LGM_fit,
                    X.pop = X_pop_sf,
                    domain = ~admin1_name + admin2_name,
                    mesh,
                    n.sample = 1000,
                    cluster.effect = F,
                    X.pop.weights = X_pop$adm2_pop_weight,
                    level = .8,
                    return.samples = T)
    bbin_LGM_res$betabinomial.spde.lgm.est.subdomain$method <- "Betabinomial GRF"
    holdout_i <- match(holdout_name,
                       bbin_LGM_res$betabinomial.spde.lgm.est.subdomain$domain)
    spde_logit_samples <-
      SUMMER::logit(bbin_LGM_res$betabinomial.spde.lgm.sample.subdomain)
    bbin_LGM_res$betabinomial.spde.lgm.est.subdomain$logit_mean <- 
      apply(spde_logit_samples, 1, mean)
    bbin_LGM_res$betabinomial.spde.lgm.est.subdomain$logit_var <- 
      apply(spde_logit_samples, 1, var)
    res_list <- 
      c(res_list,
        list(bbin_LGM_res$betabinomial.spde.lgm.est.subdomain[holdout_i, ]))
    
    holdout_res <- bind_rows(res_list)
    saveRDS(holdout_res, paste0("results/estimates/holdout_res_", holdout_name, ".rds"))
    
  } else {
    holdout_res <- readRDS(paste0("results/estimates/holdout_res_", holdout_name, ".rds"))
  }
  full_res_list <- c(full_res_list, list(holdout_res))
}

holdout_res <- bind_rows(full_res_list)

full_res_list <- list()
for (holdout_name in poly_adm2$NAME_2) {
  if (file.exists(paste0("results/estimates/holdout_res_", holdout_name, ".rds"))) {
    holdout_res <- readRDS(paste0("results/estimates/holdout_res_", holdout_name, ".rds"))
    full_res_list <- c(full_res_list, list(holdout_res))
  }
}

holdout_res <- bind_rows(full_res_list)

sample_des <- svydesign(id = ~cluster + hshold,
                        strata = ~stratum, nest = T,
                        weights = ~wt, data = svy_dat)

adm2_alm_res <- smoothArea(hiv~1, domain = ~admin2_name, 
                           design = sample_des, 
                           adj.mat = admin2_mat, 
                           transform = "logit")
holdout_res <- holdout_res |>
  left_join(adm2_alm_res$direct.est |>
              select(domain, median, var) |>
              rename(direct.est = median, 
                     direct.est.var = var),
            by = "domain") |>
  filter(!is.na(direct.est))

holdout_res <- holdout_res |>
  mutate(logit_direct_est = SUMMER::logit(direct.est),
         logit_direct_est_var = direct.est.var / direct.est^2 / (1-direct.est) ^2) |>
  mutate(logit_upper = logit_mean + qnorm(.9) * sqrt(logit_direct_est_var + logit_var),
         logit_lower = logit_mean - qnorm(.9) * sqrt(logit_direct_est_var + logit_var),
         sq_err = (median - direct.est)^2,
         abs_err = abs(median - direct.est),
         cov = logit_lower < logit_direct_est & logit_upper > logit_direct_est,
         log_cpo = dnorm(logit_direct_est, mean = logit_mean, sd = sqrt(logit_direct_est_var + logit_var), log = T),
         is = logit_upper - logit_lower + 
           2 / 0.2 * ifelse(logit_direct_est < logit_lower, logit_lower - logit_direct_est, 0) +
           2 / 0.2 * ifelse(logit_direct_est > logit_upper, logit_direct_est - logit_upper, 0))
saveRDS(holdout_res, "results/estimates/holdout_res.rds")
holdout_res <- readRDS("results/estimates/holdout_res.rds")

comp_table <- holdout_res |>
  group_by(method) |>
  summarize(mae = mean(abs_err),
            rmse = sqrt(mean(sq_err)),
            cov = mean(cov),
            log_cpo = mean(log_cpo),
            is = mean(is),
            int_length = mean(upper - lower))

present_comp_table <- comp_table |>
  filter(method %in% c("Direct",
                         "Binomial BYM2",
                         "Betabinomial BYM2",
                         "Fay-Herriot BYM2",
                         "Binomial GRF",
                         "Betabinomial GRF")) |>
  mutate(mae = round(mae * 100, 2),
         rmse = round(rmse * 100, 2),
         cov = round(cov * 100),
         log_cpo = round(log_cpo - max(log_cpo), 2),
         is = round(is, 2),
         int_length = round(int_length * 100, 2))
holdout_res |>
  filter(!stringr::str_detect(holdout_res$method, "No cov.")) |>
  ggplot(aes(x = logit_direct_est, y = log_cpo, color = method)) + 
  geom_point() + 
  facet_wrap(~method) +
  theme_bw() +
  ggtitle(label = "Log(CPO) vs logit(direct estimate)", 
          subtitle = "Large outliers for betabinomial models?")
ggsave("results/figures/Zambia_adm2_cv_log_cpo_detailed.pdf", width = 11, height = 8)
holdout_res |>
  group_by(domain) |>
  mutate(mean_log_cpo = mean(log_cpo)) |>
  ungroup() |>
  arrange(mean_log_cpo) |>
  mutate(domain = factor(domain, levels = unique(domain))) |>
  filter(!stringr::str_detect(holdout_res$method, "No cov.")) |>
  filter(domain %in% unique(domain)[1:20]) |>
  ggplot(aes(x = method, y = logit_mean, color = method)) +
  geom_point() +
  geom_hline(aes(yintercept = logit_direct_est)) + 
  geom_hline(yintercept = logit(svymean(~hiv, sample_des)[1]), 
             color = "red", linetype = "dashed") +
  geom_linerange(aes(x = method, ymin = logit_lower, ymax = logit_upper)) +
  facet_wrap(~domain) +
  scale_x_discrete(labels = NULL, breaks = NULL) +
  theme_bw() +
  ggtitle(label = "20 Admin-2 areas with largest log(CPO) values", 
          subtitle = "Red line is national mean, black line is Admin-2 direct estimate")
ggsave("results/figures/Zambia_adm2_cv_detailed.pdf", width = 11, height = 8)
holdout_res |> 
  group_by(method) |> 
  summarize(rmse = sqrt(mean((median -direct.est) ^ 2)),
            mae = mean(abs(median - direct.est)),
            cov = mean(lower < direct.est & upper > direct.est))

