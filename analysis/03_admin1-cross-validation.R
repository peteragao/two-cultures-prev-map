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

library(stringr)
geo_file <- "data/Zambia/DHS/ZMGE71FL/"
ir_file <- "data/Zambia/DHS/ZMIR71DT/ZMIR71FL.DTA"
hiv_file <- "data/Zambia/DHS/ZMAR71DT/ZMAR71FL.DTA"
gadm_abbrev <- "ZMB"
poly_path <- "data/Zambia/GADM/gadm41_ZMB_shp"

source("~/prevalence-mapping-review/analysis/functions.R")
setwd("~/prevalence-mapping-review/")
source("analysis/functions.R")
inla.setOption("num.threads", 8)
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

X_pop_sf <- st_as_sf(X_pop, coords = c("LONGNUM", "LATNUM"), crs = "EPSG:4326")

# DHS RESPONSE DATA ------------------------------------------------------------
svy_dat <- readRDS("data/Zambia/clean_DHS_data.rds")

# CROSS-VALIDATION -------------------------------------------------------------
m <- nrow(admin1_mat)
res_list <- list()
options(survey.lonely.psu="adjust")

sample_locs <- st_sample(poly_adm0, 1000)
mesh = inla.mesh.2d(loc.domain = st_coordinates(sample_locs),
                    max.edge = 0.25, offset = -0.1)
prior.range = c(3, 0.5)
prior.sigma = c(0.5, 0.5)
prior.clust  = list(prec = list(prior = "pc.prec",
                                param = c(1, 0.05)))

adm1_est_table <- read.csv("results/estimates/adm1_est.csv")


for (i in 1:10) {
  
  holdout_name <- rownames(admin1_mat)[i]
  cat("\nTake out region: ")
  cat(i, " (", holdout_name, ")", sep = "")
  cat("\n")
  print(Sys.time())
  
  holdout_dat <- svy_dat |>
    filter(admin1 != i)
  holdout_des <- svydesign(id = ~cluster + hshold,
                           strata = ~stratum, nest = T,
                           weights = ~wt, data = holdout_dat)
  #### MODELS ####
  cov_formula <- hiv ~ urban + l1pntl + access + malaria
  #### AREA LEVEL ####
  adm1_alm_res <- 
    smoothArea(hiv ~ l1pntl_adm1 + access_adm1 + malaria_adm1,
               domain = ~admin1_name, 
               X.domain = X_adm1_avg,
               design = holdout_des, 
               adj.mat = admin1_mat, 
               transform = "logit", 
               level = .8,
               return.samples = T)
  
  holdout_i <- match(holdout_name, adm1_alm_res$iid.model.est$domain)
  iid_logit_samples <-
    SUMMER::logit(adm1_alm_res$iid.model.sample)
  adm1_alm_res$iid.model.est$logit_mean <- 
    apply(iid_logit_samples, 1, mean)
  adm1_alm_res$iid.model.est$logit_var <- 
    apply(iid_logit_samples, 1, var)
  adm1_alm_res$iid.model.est$method <- "Fay-Herriot IID"
  bym2_logit_samples <-
    SUMMER::logit(adm1_alm_res$bym2.model.sample)[match(unique(adm1_est_table$domain), rownames(admin1_mat)), ]
  adm1_alm_res$bym2.model.est$logit_mean <- 
    apply(bym2_logit_samples, 1, mean)
  adm1_alm_res$bym2.model.est$logit_var <- 
    apply(bym2_logit_samples, 1, var)
  adm1_alm_res$bym2.model.est$method <- "Fay-Herriot BYM2"
  res_list <- 
    c(res_list,
      list(adm1_alm_res$iid.model.est[holdout_i, ]),
      list(adm1_alm_res$bym2.model.est[holdout_i, ]))
  
  #### AREA LEVEL NO COVARIATES ####
  adm1_alm_no_cov_res <- 
    smoothArea(hiv~1, domain = ~admin1_name, 
               design = holdout_des, 
               adj.mat = admin1_mat, 
               transform = "logit", 
               level = .8,
               return.samples = T)
  iid_logit_samples <-
    SUMMER::logit(adm1_alm_no_cov_res$iid.model.sample)
  adm1_alm_no_cov_res$iid.model.est$logit_mean <- 
    apply(iid_logit_samples, 1, mean)
  adm1_alm_no_cov_res$iid.model.est$logit_var <- 
    apply(iid_logit_samples, 1, var)
  adm1_alm_no_cov_res$iid.model.est$method <- "Fay-Herriot IID: No cov."
  bym2_logit_samples <-
    SUMMER::logit(adm1_alm_no_cov_res$bym2.model.sample)[match(unique(adm1_est_table$domain), rownames(admin1_mat)), ]
  adm1_alm_no_cov_res$bym2.model.est$logit_mean <- 
    apply(bym2_logit_samples, 1, mean)
  adm1_alm_no_cov_res$bym2.model.est$logit_var <- 
    apply(bym2_logit_samples, 1, var)
  adm1_alm_no_cov_res$bym2.model.est$method <- "Fay-Herriot BYM2: No cov."
  res_list <- 
    c(res_list,
      list(adm1_alm_no_cov_res$iid.model.est[holdout_i, ]),
      list(adm1_alm_no_cov_res$bym2.model.est[holdout_i, ]))
  
  #### BINOMIAL ####
  bin_adm1_res <-
    smoothUnit(
      cov_formula,
      domain = ~admin1_name, 
      design = holdout_des, 
      adj.mat = admin1_mat,
      X.pop = X_pop,
      family = "binomial",
      X.pop.weights = X_pop$adm1_pop_weight,
      level = .8,
      return.samples = T
    )
  bin_adm1_res$bym2.model.est$method <- "Binomial BYM2"
  holdout_i <- match(holdout_name, bin_adm1_res$bym2.model.est$domain)
  bym2_logit_samples <-
    SUMMER::logit(bin_adm1_res$bym2.model.sample)[match(unique(adm1_est_table$domain), rownames(admin1_mat)), ]
  bin_adm1_res$bym2.model.est$logit_mean <- 
    apply(bym2_logit_samples, 1, mean)
  bin_adm1_res$bym2.model.est$logit_var <- 
    apply(bym2_logit_samples, 1, var)
  res_list <- 
    c(res_list,
      list(bin_adm1_res$bym2.model.est[holdout_i, ]))
  
  #### BINOMIAL NO COVARIATES ####
  bin_no_cov_adm1_res <-
    smoothUnit(
      hiv ~ urban,
      domain = ~admin1_name, 
      design = holdout_des, 
      adj.mat = admin1_mat,
      X.pop = X_pop,
      family = "binomial",
      X.pop.weights = X_pop$adm1_pop_weight,
      level = .8,
      return.samples = T
    )
  bin_no_cov_adm1_res$bym2.model.est$method <- "Binomial BYM2: No cov."
  holdout_i <- match(holdout_name, bin_no_cov_adm1_res$bym2.model.est$domain)
  bym2_logit_samples <-
    SUMMER::logit(bin_no_cov_adm1_res$bym2.model.sample)[match(unique(adm1_est_table$domain), rownames(admin1_mat)), ]
  bin_no_cov_adm1_res$bym2.model.est$logit_mean <- 
    apply(bym2_logit_samples, 1, mean)
  bin_no_cov_adm1_res$bym2.model.est$logit_var <- 
    apply(bym2_logit_samples, 1, var)
  res_list <- 
    c(res_list,
      list(bin_no_cov_adm1_res$bym2.model.est[holdout_i, ]))
  
  
  #### BETABINOMIAL ####
  bbin_adm1_res <-
    fitBetabinomialBYM2(
      cov_formula,
      domain = ~admin1_name, 
      design = holdout_des, 
      adj.mat = admin1_mat,
      X.pop = X_pop,
      X.pop.weights = X_pop$adm1_pop_weight,
      level = .8,
      return.samples = T
    )
  bbin_adm1_res$betabinomial.bym2.model.est$method <- "Betabinomial BYM2"
  holdout_i <- match(holdout_name, bbin_adm1_res$betabinomial.bym2.model.est$domain)
  bym2_logit_samples <-
    SUMMER::logit(bbin_adm1_res$betabinomial.bym2.model.sample)[match(unique(adm1_est_table$domain), rownames(admin1_mat)), ]
  bbin_adm1_res$betabinomial.bym2.model.est$logit_mean <- 
    apply(bym2_logit_samples, 1, mean)
  bbin_adm1_res$betabinomial.bym2.model.est$logit_var <- 
    apply(bym2_logit_samples, 1, var)
  res_list <- 
    c(res_list,
      list(bbin_adm1_res$betabinomial.bym2.model.est[holdout_i, ]))
  #### BETABINOMIAL NO COVARIATES ####
  bbin_no_cov_adm1_res <-
    fitBetabinomialBYM2(
      hiv ~ urban,
      domain = ~admin1_name, 
      design = holdout_des, 
      adj.mat = admin1_mat,
      X.pop = X_pop,
      X.pop.weights = X_pop$adm1_pop_weight,
      level = .8,
      return.samples = T
    )
  bbin_no_cov_adm1_res$betabinomial.bym2.model.est$method <-
    "Betabinomial BYM2: No cov."
  holdout_i <- match(holdout_name, bbin_no_cov_adm1_res$betabinomial.bym2.model.est$domain)
  bym2_logit_samples <-
    SUMMER::logit(bbin_no_cov_adm1_res$betabinomial.bym2.model.sample)[match(unique(adm1_est_table$domain), rownames(admin1_mat)), ]
  bbin_no_cov_adm1_res$betabinomial.bym2.model.est$logit_mean <- 
    apply(bym2_logit_samples, 1, mean)
  bbin_no_cov_adm1_res$betabinomial.bym2.model.est$logit_var <- 
    apply(bym2_logit_samples, 1, var)
  res_list <- 
    c(res_list,
      list(bbin_no_cov_adm1_res$betabinomial.bym2.model.est[holdout_i, ]))
  
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
                  X.pop.weights = X_pop$adm1_pop_weight,
                  level = .8,
                  return.samples = T)
  bin_LGM_res$binomial.spde.lgm.est$method <- "Binomial GRF"
  holdout_i <- match(holdout_name,  bin_LGM_res$binomial.spde.lgm.est$domain)
  
  spde_logit_samples <-
    SUMMER::logit(bin_LGM_res$binomial.spde.lgm.sample)
  bin_LGM_res$binomial.spde.lgm.est$logit_mean <- 
    apply(spde_logit_samples, 1, mean)
  bin_LGM_res$binomial.spde.lgm.est$logit_var <- 
    apply(spde_logit_samples, 1, var)
  
  res_list <- 
    c(res_list,
      list(bin_LGM_res$binomial.spde.lgm.est[holdout_i, ]))
  
  
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
                  X.pop.weights = X_pop$adm1_pop_weight,
                  level = .8,
                  return.samples = T)
  bbin_LGM_res$betabinomial.spde.lgm.est$method <- "Betabinomial GRF"
  holdout_i <- match(holdout_name,
                     bbin_LGM_res$betabinomial.spde.lgm.est$domain)
  spde_logit_samples <-
    SUMMER::logit(bbin_LGM_res$betabinomial.spde.lgm.sample)
  bbin_LGM_res$betabinomial.spde.lgm.est$logit_mean <- 
    apply(spde_logit_samples, 1, mean)
  bbin_LGM_res$betabinomial.spde.lgm.est$logit_var <- 
    apply(spde_logit_samples, 1, var)
  res_list <- 
    c(res_list,
      list(bbin_LGM_res$betabinomial.spde.lgm.est[holdout_i, ]))
  #### BETABINOMIAL SPDE NO COV ####
  bbin_LGM_no_cov_fit <- 
    fitContLGM(formula = hiv ~ urban,
               family = "betabinomial",
               cluster = ~cluster,
               cluster.effect = F,
               data = holdout_dat,
               mesh = mesh,
               pc.prior.range = prior.range,
               pc.prior.sigma = prior.sigma,
               pc.prior.clust = prior.clust)
  bbin_LGM_no_cov_res <- 
    smoothContLGM(bbin_LGM_no_cov_fit,
                  X.pop = X_pop_sf,
                  domain = ~admin1_name + admin2_name,
                  mesh,
                  n.sample = 1000,
                  cluster.effect = F,
                  X.pop.weights = X_pop$adm1_pop_weight,
                  level = .8,
                  return.samples = T)
  bbin_LGM_no_cov_res$betabinomial.spde.lgm.est$method <- "Betabinomial GRF: No cov."
  holdout_i <- match(holdout_name,
                     bbin_LGM_no_cov_res$betabinomial.spde.lgm.est$domain)
  spde_logit_samples <-
    SUMMER::logit(bbin_LGM_no_cov_res$betabinomial.spde.lgm.sample)
  bbin_LGM_no_cov_res$betabinomial.spde.lgm.est$logit_mean <- 
    apply(spde_logit_samples, 1, mean)
  bbin_LGM_no_cov_res$betabinomial.spde.lgm.est$logit_var <- 
    apply(spde_logit_samples, 1, var)
  res_list <- 
    c(res_list,
      list(bbin_LGM_no_cov_res$betabinomial.spde.lgm.est[holdout_i, ]))
 
}

holdout_res <- bind_rows(res_list)
saveRDS(holdout_res, "results/estimates/adm1_holdout_res.rds")


holdout_res <- readRDS("results/estimates/adm1_holdout_res.rds")
res_list <- list(holdout_res)
for (i in 1:10) {
  
  holdout_name <- rownames(admin1_mat)[i]
  cat("\nTake out region: ")
  cat(i, " (", holdout_name, ")", sep = "")
  cat("\n")
  print(Sys.time())
  
  holdout_dat <- svy_dat |>
    filter(admin1 != i)
  holdout_des <- svydesign(id = ~cluster + hshold,
                           strata = ~stratum, nest = T,
                           weights = ~wt, data = holdout_dat)
  #### BETABINOMIAL SPDE NO COV ####
  print("GRF MODELS")
  bbin_LGM_no_cov_fit <- 
    fitContLGM(formula = hiv ~ urban,
               family = "betabinomial",
               cluster = ~cluster,
               cluster.effect = F,
               data = holdout_dat,
               mesh = mesh,
               pc.prior.range = prior.range,
               pc.prior.sigma = prior.sigma,
               pc.prior.clust = prior.clust)
  bbin_LGM_no_cov_res <- 
    smoothContLGM(bbin_LGM_no_cov_fit,
                  X.pop = X_pop_sf,
                  domain = ~admin1_name + admin2_name,
                  mesh,
                  n.sample = 1000,
                  cluster.effect = F,
                  X.pop.weights = X_pop$adm1_pop_weight,
                  level = .8,
                  return.samples = T)
  bbin_LGM_no_cov_res$betabinomial.spde.lgm.est$method <- "Betabinomial GRF: No cov."
  holdout_i <- match(holdout_name,
                     bbin_LGM_no_cov_res$betabinomial.spde.lgm.est$domain)
  spde_logit_samples <-
    SUMMER::logit(bbin_LGM_no_cov_res$betabinomial.spde.lgm.sample)
  bbin_LGM_no_cov_res$betabinomial.spde.lgm.est$logit_mean <- 
    apply(spde_logit_samples, 1, mean)
  bbin_LGM_no_cov_res$betabinomial.spde.lgm.est$logit_var <- 
    apply(spde_logit_samples, 1, var)
  res_list <- 
    c(res_list,
      list(bbin_LGM_no_cov_res$betabinomial.spde.lgm.est[holdout_i, ]))
  
}
