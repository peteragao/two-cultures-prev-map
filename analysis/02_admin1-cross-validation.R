# 02_admin1-cross-validation.R


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

geo_file <- "data/Zambia/DHS/ZMGE71FL/"
ir_file <- "data/Zambia/DHS/ZMIR71DT/ZMIR71FL.DTA"
hiv_file <- "data/Zambia/DHS/ZMAR71DT/ZMAR71FL.DTA"

source("~/prevalence-mapping-review/analysis/functions.R")
setwd("~/prevalence-mapping-review/")
#INLA:::inla.binary.install()
inla.setOption("num.threads", 16)
# GADM POLYGON DATA ------------------------------------------------------------
# GADM version 3.6 is used to reflect Admin-1 and Admin-2 boundaries
gadm_abbrev <- "ZMB"
poly_path <- "data/Zambia/GADM/gadm36_ZMB_shp"
poly_layer_adm0 <- paste('gadm36', gadm_abbrev,
                         '0', sep = "_")
poly_layer_adm1 <- paste('gadm36', gadm_abbrev,
                         '1', sep = "_")
poly_layer_adm2 <- paste('gadm36', gadm_abbrev,
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

# transform
X_pop$l1pntl <- log1p(X_pop$ntl)
X_pop$l1paccess <- log1p(X_pop$access)
X_pop$urban <- 1 * (X_pop$urban)

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
  select(admin1_name, l1pntl_adm1, l1paccess_adm1, malaria_adm1) |>
  group_by(admin1_name) |>
  slice(which.min(l1pntl_adm1)) |>
  ungroup() |>
  st_set_geometry(NULL)
X_adm2_avg <- X_pop |>
  select(admin2_name, l1pntl_adm2, l1paccess_adm2, malaria_adm2) |>
  group_by(admin2_name) |>
  slice(which.min(l1pntl_adm2)) |>
  ungroup() |>
  st_set_geometry(NULL)

# DHS RESPONSE DATA ------------------------------------------------------------
svy_dat <- readRDS("data/Zambia/clean_DHS_data.rds")
svy_dat <- svy_dat |>
  left_join(X_adm1_avg, by = "admin1_name") |>
  left_join(X_adm2_avg, by = "admin2_name")
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
  #### AREA LEVEL ####
  adm1_alm_res <- 
    smoothArea(hiv ~ l1pntl_adm1 + l1paccess_adm1 + malaria_adm1,
               domain = ~admin1_name, 
               X.domain = X_adm1_avg,
               design = holdout_des, 
               adj.mat = admin1_mat, 
               transform = "logit", 
               level = .5)

  holdout_i <- match(holdout_name, adm1_alm_res$iid.model.est$domain)
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
               level = .5)
  res_list <- 
    c(res_list,
      list(adm1_alm_no_cov_res$iid.model.est[holdout_i, ] |>
             mutate(method = paste0(method, ": No cov."))),
      list(adm1_alm_no_cov_res$bym2.model.est[holdout_i, ] |>
             mutate(method = paste0(method, ": No cov."))))

  #### BINOMIAL ####
  bin_adm1_res <-
    smoothUnit(
      hiv ~ urban + l1pntl*urban + l1paccess*urban + malaria*urban,
      domain = ~admin1_name, 
      design = holdout_des, 
      adj.mat = admin1_mat,
      X.pop = X_pop,
      family = "binomial",
      X.pop.weights = X_pop$adm1_pop_weight,
      level = .5
    )
  bin_adm1_res$bym2.model.est$method <- "Binomial BYM2"
  holdout_i <- match(holdout_name, bin_adm1_res$bym2.model.est$domain)
  
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
      level = .5
    )
  bin_no_cov_adm1_res$bym2.model.est$method <- "Binomial BYM2: No cov."
  holdout_i <- match(holdout_name, bin_no_cov_adm1_res$bym2.model.est$domain)
  res_list <- 
    c(res_list,
      list(bin_no_cov_adm1_res$bym2.model.est[holdout_i, ]))
  
  #### LONOBINOMIAL ####
  lono_adm1_res <- 
    fitLonobinomialBYM2(
      hiv ~ urban +  l1pntl*urban + l1paccess*urban + malaria*urban,
      domain = ~admin1_name, 
      design = holdout_des, 
      adj.mat = admin1_mat,
      X.pop = X_pop,
      X.pop.weights = X_pop$adm1_pop_weight,
      level = .5
    )
  lono_adm1_res$lonobinomial.bym2.model.est$method <- "Lonobinomial BYM2"
  holdout_i <- match(holdout_name, lono_adm1_res$lonobinomial.bym2.model.est$domain)
  res_list <- 
    c(res_list,
      list(lono_adm1_res$lonobinomial.bym2.model.est[holdout_i, ]))
  
  #### LONOBINOMIAL NO COVARIATES ####
  lono_no_cov_adm1_res <- 
    fitLonobinomialBYM2(
      hiv ~ urban,
      domain = ~admin1_name, 
      design = holdout_des, 
      adj.mat = admin1_mat,
      X.pop = X_pop,
      X.pop.weights = X_pop$adm1_pop_weight,
      level = .5
    )
  lono_no_cov_adm1_res$lonobinomial.bym2.model.est$method <-
    "Lonobinomial BYM2: No cov."
  holdout_i <- match(holdout_name, lono_no_cov_adm1_res$lonobinomial.bym2.model.est$domain)
  res_list <- 
    c(res_list,
      list(lono_no_cov_adm1_res$lonobinomial.bym2.model.est[holdout_i, ]))
  
  #### BETABINOMIAL ####
  bbin_adm1_res <-
    fitBetabinomialBYM2(
      hiv ~ urban + l1pntl*urban + l1paccess*urban + malaria*urban,
      domain = ~admin1_name, 
      design = holdout_des, 
      adj.mat = admin1_mat,
      X.pop = X_pop,
      X.pop.weights = X_pop$adm1_pop_weight,
      level = .5
    )
  bbin_adm1_res$betabinomial.bym2.model.est$method <- "Betabinomial BYM2"
  holdout_i <- match(holdout_name, bbin_adm1_res$betabinomial.bym2.model.est$domain)
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
      level = .5
    )
  bbin_no_cov_adm1_res$betabinomial.bym2.model.est$method <-
    "Betabinomial BYM2: No cov."
  holdout_i <- match(holdout_name, bbin_no_cov_adm1_res$betabinomial.bym2.model.est$domain)
  res_list <- 
    c(res_list,
      list(bbin_no_cov_adm1_res$betabinomial.bym2.model.est[holdout_i, ]))
  
  #### BINOMIAL SPDE ####
  bin_LGM_fit <- 
    fitContLGM(formula = hiv ~ urban + l1pntl*urban + l1paccess*urban + malaria*urban,
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
                  X.pop = X_pop,
                  domain = ~admin1_name + admin2_name,
                  mesh,
                  n.sample = 1000,
                  cluster.effect = F,
                  X.pop.weights = X_pop$adm1_pop_weight,
                  level = 0.5)
  bin_LGM_res$binomial.spde.lgm.est$method <- "Binomial SPDE LGM"
  holdout_i <- match(holdout_name,  bin_LGM_res$binomial.spde.lgm.est$domain)
  res_list <- 
    c(res_list,
      list(bin_LGM_res$binomial.spde.lgm.est[holdout_i, ]))
  #### BINOMIAL SPDE NO COVARIATES ####
  bin_LGM_no_cov_fit <- 
    fitContLGM(formula = hiv ~ urban,
               family = "binomial",
               cluster = ~cluster,
               cluster.effect = F,
               data = holdout_dat,
               mesh = mesh,
               pc.prior.range = prior.range,
               pc.prior.sigma = prior.sigma,
               pc.prior.clust = prior.clust)
  bin_LGM_no_cov_res <- 
    smoothContLGM(bin_LGM_no_cov_fit,
                  X.pop = X_pop,
                  domain = ~admin1_name + admin2_name,
                  mesh,
                  n.sample = 1000,
                  cluster.effect = F,
                  X.pop.weights = X_pop$adm1_pop_weight,
                  level = 0.5)
  bin_LGM_no_cov_res$binomial.spde.lgm.est$method <- 
    "Binomial SPDE LGM: No cov."
  holdout_i <- match(holdout_name, 
                     bin_LGM_no_cov_res$binomial.spde.lgm.est$domain)
  res_list <- 
    c(res_list,
      list(bin_LGM_no_cov_res$binomial.spde.lgm.est[holdout_i, ]))
  
  #### LONOBINOMIAL SPDE ####
  lono_LGM_fit <- 
    fitContLGM(formula = hiv ~ urban + l1pntl*urban + l1paccess*urban + malaria*urban,
               family = "binomial",
               cluster = ~cluster,
               cluster.effect = T,
               data = holdout_dat,
               mesh = mesh,
               pc.prior.range = prior.range,
               pc.prior.sigma = prior.sigma,
               pc.prior.clust = prior.clust)
  lono_LGM_res <- 
    smoothContLGM(lono_LGM_fit,
                  X.pop = X_pop,
                  domain = ~admin1_name + admin2_name,
                  mesh,
                  n.sample = 1000,
                  cluster.effect = T,
                  X.pop.weights = X_pop$adm1_pop_weight,
                  level = 0.5)
  lono_LGM_res$binomial.spde.lgm.est$method <- "Lonobinomial SPDE LGM"
  holdout_i <- match(holdout_name, lono_LGM_res$binomial.spde.lgm.est$domain)
  res_list <- 
    c(res_list,
      list(lono_LGM_res$binomial.spde.lgm.est[holdout_i, ]))
  
  #### LONOBINOMIAL SPDE NO COVARIATES ####
  lono_LGM_no_cov_fit <- 
    fitContLGM(formula = hiv ~ urban,
               family = "binomial",
               cluster = ~cluster,
               cluster.effect = T,
               data = holdout_dat,
               mesh = mesh,
               pc.prior.range = prior.range,
               pc.prior.sigma = prior.sigma,
               pc.prior.clust = prior.clust)
  lono_LGM_no_cov_res <- 
    smoothContLGM(lono_LGM_no_cov_fit,
                  X.pop = X_pop,
                  domain = ~admin1_name + admin2_name,
                  mesh,
                  n.sample = 1000,
                  cluster.effect = T,
                  X.pop.weights = X_pop$adm1_pop_weight,
                  level = 0.5)
  lono_LGM_no_cov_res$binomial.spde.lgm.est$method <-
    "Lonobinomial SPDE LGM: No cov."
  holdout_i <- match(holdout_name,
                     lono_LGM_no_cov_res$binomial.spde.lgm.est$domain)
  res_list <- 
    c(res_list,
      list(lono_LGM_no_cov_res$binomial.spde.lgm.est[holdout_i, ]))
  
  #### BETABINOMIAL SPDE ####
  bbin_LGM_fit <- 
    fitContLGM(formula = hiv ~ urban + l1pntl*urban + l1paccess*urban + malaria*urban,
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
                  X.pop = X_pop,
                  domain = ~admin1_name + admin2_name,
                  mesh,
                  n.sample = 1000,
                  cluster.effect = F,
                  X.pop.weights = X_pop$adm1_pop_weight,
                  level = 0.5)
  bbin_LGM_res$betabinomial.spde.lgm.est$method <- "Betabinomial SPDE LGM"
  holdout_i <- match(holdout_name,
                     bbin_LGM_res$betabinomial.spde.lgm.est$domain)
  res_list <- 
    c(res_list,
      list(bbin_LGM_res$betabinomial.spde.lgm.est[holdout_i, ]))
  #### BETABINOMIAL SPDE NO COVARIATES ####
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
                  X.pop = X_pop,
                  domain = ~admin1_name + admin2_name,
                  mesh,
                  n.sample = 1000,
                  cluster.effect = F,
                  X.pop.weights = X_pop$adm1_pop_weight,
                  level = 0.5)
  bbin_LGM_no_cov_res$betabinomial.spde.lgm.est$method <- "Betabinomial SPDE LGM: No cov."
  holdout_i <- match(holdout_name,
                     bbin_LGM_no_cov_res$betabinomial.spde.lgm.est$domain)
  res_list <- 
    c(res_list,
      list(bbin_LGM_no_cov_res$betabinomial.spde.lgm.est[holdout_i, ]))
  
  #### LONOBINOMIAL BYM2 FRAME ####
  lono_simframe_adm1_res <- 
    fitLonobinomialBYM2(
      hiv ~ urban +  l1pntl*urban + l1paccess*urban + malaria*urban,
      domain = ~admin1_name, 
      design = holdout_des, 
      adj.mat = admin1_mat,
      X.pop = X_sim_frame,
      X.pop.weights = X_sim_frame$adm1_pop_weight,
      level = .5
    )
  lono_simframe_adm1_res$lonobinomial.bym2.model.est$method <-
    "Lonobinomial BYM2 Sim. Frame"
  holdout_i <- match(holdout_name, 
                     lono_simframe_adm1_res$lonobinomial.bym2.model.est$domain)
  res_list <- 
    c(res_list,
      list(lono_simframe_adm1_res$lonobinomial.bym2.model.est[holdout_i, ]))
  #### BETABINOMIAL BYM2 FRAME ####
  bbin_simframe_adm1_res <-
    fitBetabinomialBYM2(
      hiv ~ urban + l1pntl*urban + l1paccess*urban + malaria*urban,
      domain = ~admin1_name, 
      design = holdout_des, 
      adj.mat = admin1_mat,
      X.pop = X_sim_frame,
      X.pop.weights = X_sim_frame$adm1_pop_weight,
      level = .5
    )
  bbin_simframe_adm1_res$betabinomial.bym2.model.est$method <-
    "Betabinomial BYM2 Sim. Frame"
  holdout_i <- match(holdout_name, 
                     bbin_simframe_adm1_res$betabinomial.bym2.model.est$domain)
  res_list <- 
    c(res_list,
      list(bbin_simframe_adm1_res$betabinomial.bym2.model.est[holdout_i, ]))
  #### LONOBINOMIAL SPDE FRAME ####
  lono_LGM_simframe_res <- 
    smoothContLGM(lono_LGM_fit,
                  X.pop = X_sim_frame,
                  domain = ~admin1_name + admin2_name,
                  mesh,
                  n.sample = 1000,
                  cluster.effect = T,
                  X.pop.weights = X_sim_frame$adm1_pop_weight,
                  level = 0.5)
  lono_LGM_simframe_res$binomial.spde.lgm.est$method <-
    "Lonobinomial SPDE LGM Sim. Frame"
  holdout_i <- match(holdout_name, 
                     lono_LGM_simframe_res$binomial.spde.lgm.est$domain)
  res_list <- 
    c(res_list,
      list(lono_LGM_simframe_res$binomial.spde.lgm.est[holdout_i, ]))
  #### BETABINOMIAL SPDE FRAME ####
  bbin_LGM_simframe_res <- 
    smoothContLGM(bbin_LGM_fit,
                  X.pop = X_sim_frame,
                  domain = ~admin1_name + admin2_name,
                  mesh,
                  n.sample = 1000,
                  cluster.effect = F,
                  X.pop.weights = X_sim_frame$adm1_pop_weight,
                  level = 0.5)
  bbin_LGM_simframe_res$betabinomial.spde.lgm.est$method <-
    "Betabinomial SPDE LGM Sim. Frame"
  holdout_i <- match(holdout_name,
                     bbin_LGM_simframe_res$betabinomial.spde.lgm.est$domain)
  res_list <- 
    c(res_list,
      list(bbin_LGM_simframe_res$betabinomial.spde.lgm.est[holdout_i, ]))
  
}

holdout_res <- bind_rows(res_list)

sample_des <- svydesign(id = ~cluster + hshold,
                        strata = ~stratum, nest = T,
                        weights = ~wt, data = svy_dat)

adm1_alm_res <- smoothArea(hiv~1, domain = ~admin1_name, 
                           design = sample_des, 
                           adj.mat = admin1_mat, 
                           transform = "logit")
holdout_res <- holdout_res |>
  left_join(adm1_alm_res$direct.est |>
              select(domain, median, var) |>
              rename(direct.est = median, 
                     direct.est.var = var),
            by = "domain") 

saveRDS(holdout_res, "results/estimates/holdout_res_adm1.rds")

holdout_res |> 
  group_by(method) |> 
  summarize(rmse = sqrt(mean((median -direct.est) ^ 2)),
            mae = mean(abs(median - direct.est)),
            cov = mean(lower < direct.est & upper > direct.est))






