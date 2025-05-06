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

geo_file <- "data/Zambia/DHS/ZMGE71FL/"
ir_file <- "data/Zambia/DHS/ZMIR71DT/ZMIR71FL.DTA"
hiv_file <- "data/Zambia/DHS/ZMAR71DT/ZMAR71FL.DTA"
#source("analysis/functions.R")
source("analysis/functions.R")
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

# ID HOLDOUT AREAS -------------------------------------------------------------
# not all admin-2 areas have reliable direct estimates
# keep those with direct estimates with CV < 20% and variance > 0
adm2_est_table <- read.csv("results/estimates/adm2_est.csv")
selected_methods <- 
  c("Direct",
    "Area level BYM2",
    "Binomial BYM2",
    "Betabinomial BYM2",
    "Binomial SPDE LGM",
    "Betabinomial SPDE LGM")
reordered_methods <- 
  c("Direct",
    "Binomial BYM2",
    "Betabinomial BYM2",
    "Area level BYM2",
    "Binomial SPDE LGM",
    "Betabinomial SPDE LGM")
adm2_est_table <- adm2_est_table |>
  filter(method %in% selected_methods) |>
  mutate(method = factor(method, levels = selected_methods))
plot_sf <- poly_adm2 |>
  left_join(adm2_est_table, by = "domain")
plot_sf <- plot_sf |>
  mutate(CV = sqrt(var) / median * 100)  |>
  mutate(method = factor(method, levels = reordered_methods))
# ggplot(plot_sf) + 
#   geom_sf(aes(fill = CV < 20 & var > 1e-10)) +
#   facet_wrap(~method) +
#   theme_minimal()  + my_theme 


  

holdout_direct_est <- adm2_est_table |>
  dplyr::filter(method == "Direct") |> 
  dplyr::filter(!is.na(mean)) |>
  mutate(CV = sqrt(var) / median * 100) |>
  filter((domain %in% c("Chama", "Chipata", "Lunga", "Mpika", "Mpulungu")) & var > 1e-10)
select_levels <-
  holdout_direct_est$domain[order(holdout_direct_est$mean, decreasing = T)] 
# top 5 and bottom 5 
int_table <- adm2_est_table |>
  filter(domain %in% select_levels)
int_table$domain <- factor(int_table$domain, 
                           levels = select_levels)
int_table <- int_table |> dplyr::arrange(domain, method)
# ggplot2::ggplot(int_table, ggplot2::aes(x = domain, y = median, color = method)) +
#   ggplot2::geom_point(position = ggplot2::position_dodge(width = .6)) + 
#   ggplot2::geom_linerange(ggplot2::aes(x = domain, ymin = lower, ymax = upper), 
#                           position = ggplot2::position_dodge(width = .6)) + 
#   ggplot2::scale_color_discrete(name = "Method") + 
#   ggplot2::ylab("Small Area Estimate") + ggplot2::xlab("Domain") +
#   ggplot2::theme_bw() +
#   ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1),
#                  legend.position="bottom")

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
# easpa <- readr::read_csv("data/Zambia/DHS/zmb_2018_easpa.csv") |>
#   rename(admin1_name = GADMarea) |>
#   select(admin1_name, eaUrb, eaRur) |>
#   pivot_longer(cols = c(eaUrb, eaRur), values_to = "n_ea") |>
#   mutate(urban = 1 * (name == "eaUrb"))
# easize <- readr::read_csv("data/Zambia/DHS/zmb_2018_easpa.csv") |>
#   rename(admin1_name = GADMarea) |>
#   select(admin1_name, avgSizeUrb, avgSizeRur) |>
#   pivot_longer(cols = c(avgSizeUrb, avgSizeRur), values_to = "size") |>
#   mutate(urban = 1 * (name == "avgSizeUrb"))
# X_sim_frame <- simulateFrame(X_pop, easpa, easize)


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
      smoothArea(hiv ~ l1pntl_adm2 + l1paccess_adm2 + malaria_adm2,
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
        hiv ~ urban + l1pntl*urban + l1paccess*urban + malaria*urban,
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
    
    # #### LONOBINOMIAL ####
    # lono_adm2_res <- 
    #   fitLonobinomialBYM2(
    #     hiv ~ urban +  l1pntl*urban + l1paccess*urban + malaria*urban,
    #     domain = ~admin2_name, 
    #     design = holdout_des, 
    #     adj.mat = admin2_mat,
    #     X.pop = X_pop,
    #     X.pop.weights = X_pop$adm2_pop_weight,
    #     level = .8
    #   )
    # lono_adm2_res$lonobinomial.bym2.model.est$method <- "Lonobinomial BYM2"
    # holdout_i <- match(holdout_name, lono_adm2_res$lonobinomial.bym2.model.est$domain)
    # res_list <- 
    #   c(res_list,
    #     list(lono_adm2_res$lonobinomial.bym2.model.est[holdout_i, ]))
    # 
    # #### LONOBINOMIAL NO COVARIATES ####
    # lono_no_cov_adm2_res <- 
    #   fitLonobinomialBYM2(
    #     hiv ~ urban,
    #     domain = ~admin2_name, 
    #     design = holdout_des, 
    #     adj.mat = admin2_mat,
    #     X.pop = X_pop,
    #     X.pop.weights = X_pop$adm2_pop_weight,
    #     level = .8
    #   )
    # lono_no_cov_adm2_res$lonobinomial.bym2.model.est$method <-
    #   "Lonobinomial BYM2: No cov."
    # holdout_i <- match(holdout_name, lono_no_cov_adm2_res$lonobinomial.bym2.model.est$domain)
    # res_list <- 
    #   c(res_list,
    #     list(lono_no_cov_adm2_res$lonobinomial.bym2.model.est[holdout_i, ]))
    # 
    #### BETABINOMIAL ####
    bbin_adm2_res <-
      fitBetabinomialBYM2(
        hiv ~ urban + l1pntl*urban + l1paccess*urban + malaria*urban,
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
                    X.pop.weights = X_pop$adm2_pop_weight,
                    level = .8,
                    return.samples = T)
    bin_LGM_no_cov_res$binomial.spde.lgm.est.subdomain$method <- 
      "Binomial GRF: No cov."
    holdout_i <- match(holdout_name, 
                       bin_LGM_no_cov_res$binomial.spde.lgm.est.subdomain$domain)
    
    spde_logit_samples <-
      SUMMER::logit(bin_LGM_no_cov_res$binomial.spde.lgm.sample.subdomain)
    bin_LGM_no_cov_res$binomial.spde.lgm.est.subdomain$logit_mean <- 
      apply(spde_logit_samples, 1, mean)
    bin_LGM_no_cov_res$binomial.spde.lgm.est.subdomain$logit_var <- 
      apply(spde_logit_samples, 1, var)
    res_list <- 
      c(res_list,
        list(bin_LGM_no_cov_res$binomial.spde.lgm.est.subdomain[holdout_i, ]))
    
    # #### LONOBINOMIAL SPDE ####
    # lono_LGM_fit <- 
    #   fitContLGM(formula = hiv ~ urban + l1pntl*urban + l1paccess*urban + malaria*urban,
    #              family = "binomial",
    #              cluster = ~cluster,
    #              cluster.effect = T,
    #              data = holdout_dat,
    #              mesh = mesh,
    #              pc.prior.range = prior.range,
    #              pc.prior.sigma = prior.sigma,
    #              pc.prior.clust = prior.clust)
    # lono_LGM_res <- 
    #   smoothContLGM(lono_LGM_fit,
    #                 X.pop = X_pop,
    #                 domain = ~admin1_name + admin2_name,
    #                 mesh,
    #                 n.sample = 1000,
    #                 cluster.effect = T,
    #                 X.pop.weights = X_pop$adm2_pop_weight,
    #                 level = .8)
    # lono_LGM_res$binomial.spde.lgm.est.subdomain$method <- "Lonobinomial SPDE LGM"
    # holdout_i <- match(holdout_name, lono_LGM_res$binomial.spde.lgm.est.subdomain$domain)
    # res_list <- 
    #   c(res_list,
    #     list(lono_LGM_res$binomial.spde.lgm.est.subdomain[holdout_i, ]))
    # 
    # #### LONOBINOMIAL SPDE NO COVARIATES ####
    # lono_LGM_no_cov_fit <- 
    #   fitContLGM(formula = hiv ~ urban,
    #              family = "binomial",
    #              cluster = ~cluster,
    #              cluster.effect = T,
    #              data = holdout_dat,
    #              mesh = mesh,
    #              pc.prior.range = prior.range,
    #              pc.prior.sigma = prior.sigma,
    #              pc.prior.clust = prior.clust)
    # lono_LGM_no_cov_res <- 
    #   smoothContLGM(lono_LGM_no_cov_fit,
    #                 X.pop = X_pop,
    #                 domain = ~admin1_name + admin2_name,
    #                 mesh,
    #                 n.sample = 1000,
    #                 cluster.effect = T,
    #                 X.pop.weights = X_pop$adm2_pop_weight,
    #                 level = .8)
    # lono_LGM_no_cov_res$binomial.spde.lgm.est.subdomain$method <-
    #   "Lonobinomial SPDE LGM: No cov."
    # holdout_i <- match(holdout_name,
    #                    lono_LGM_no_cov_res$binomial.spde.lgm.est.subdomain$domain)
    # res_list <- 
    #   c(res_list,
    #     list(lono_LGM_no_cov_res$binomial.spde.lgm.est.subdomain[holdout_i, ]))
    # 
    #### BETABINOMIAL SPDE ####
    print("GRF MODELS")
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
                    X.pop.weights = X_pop$adm2_pop_weight,
                    level = .8,
                    return.samples = T)
    bbin_LGM_no_cov_res$betabinomial.spde.lgm.est.subdomain$method <- "Betabinomial GRF: No cov."
    holdout_i <- match(holdout_name,
                       bbin_LGM_no_cov_res$betabinomial.spde.lgm.est.subdomain$domain)
    spde_logit_samples <-
      SUMMER::logit(bbin_LGM_no_cov_res$betabinomial.spde.lgm.sample.subdomain)
    bbin_LGM_no_cov_res$betabinomial.spde.lgm.est.subdomain$logit_mean <- 
      apply(spde_logit_samples, 1, mean)
    bbin_LGM_no_cov_res$betabinomial.spde.lgm.est.subdomain$logit_var <- 
      apply(spde_logit_samples, 1, var)
    res_list <- 
      c(res_list,
        list(bbin_LGM_no_cov_res$betabinomial.spde.lgm.est.subdomain[holdout_i, ]))
    
    # #### LONOBINOMIAL BYM2 FRAME ####
    # lono_simframe_adm2_res <- 
    #   fitLonobinomialBYM2(
    #     hiv ~ urban +  l1pntl*urban + l1paccess*urban + malaria*urban,
    #     domain = ~admin2_name, 
    #     design = holdout_des, 
    #     adj.mat = admin2_mat,
    #     X.pop = X_sim_frame,
    #     X.pop.weights = X_sim_frame$adm2_pop_weight,
    #     level = .8
    #   )
    # lono_simframe_adm2_res$lonobinomial.bym2.model.est$method <-
    #   "Lonobinomial BYM2 Sim. Frame"
    # holdout_i <- match(holdout_name, 
    #                    lono_simframe_adm2_res$lonobinomial.bym2.model.est$domain)
    # res_list <- 
    #   c(res_list,
    #     list(lono_simframe_adm2_res$lonobinomial.bym2.model.est[holdout_i, ]))
    # #### BETABINOMIAL BYM2 FRAME ####
    # bbin_simframe_adm2_res <-
    #   fitBetabinomialBYM2(
    #     hiv ~ urban + l1pntl*urban + l1paccess*urban + malaria*urban,
    #     domain = ~admin2_name, 
    #     design = holdout_des, 
    #     adj.mat = admin2_mat,
    #     X.pop = X_sim_frame,
    #     X.pop.weights = X_sim_frame$adm2_pop_weight,
    #     level = .8
    #   )
    # bbin_simframe_adm2_res$betabinomial.bym2.model.est$method <-
    #   "Betabinomial BYM2 Sim. Frame"
    # holdout_i <- match(holdout_name, 
    #                    bbin_simframe_adm2_res$betabinomial.bym2.model.est$domain)
    # res_list <- 
    #   c(res_list,
    #     list(bbin_simframe_adm2_res$betabinomial.bym2.model.est[holdout_i, ]))
    # #### LONOBINOMIAL SPDE FRAME ####
    # lono_LGM_simframe_res <- 
    #   smoothContLGM(lono_LGM_fit,
    #                 X.pop = X_sim_frame,
    #                 domain = ~admin1_name + admin2_name,
    #                 mesh,
    #                 n.sample = 1000,
    #                 cluster.effect = T,
    #                 X.pop.weights = X_sim_frame$adm2_pop_weight,
    #                 level = .8)
    # lono_LGM_simframe_res$binomial.spde.lgm.est.subdomain$method <-
    #   "Lonobinomial SPDE LGM Sim. Frame"
    # holdout_i <- match(holdout_name, 
    #                    lono_LGM_simframe_res$binomial.spde.lgm.est.subdomain$domain)
    # res_list <- 
    #   c(res_list,
    #     list(lono_LGM_simframe_res$binomial.spde.lgm.est.subdomain[holdout_i, ]))
    # 
    # 
    # #### BETABINOMIAL SPDE FRAME ####
    # bbin_LGM_simframe_res <- 
    #   smoothContLGM(bbin_LGM_fit,
    #                 X.pop = X_sim_frame,
    #                 domain = ~admin1_name + admin2_name,
    #                 mesh,
    #                 n.sample = 1000,
    #                 cluster.effect = F,
    #                 X.pop.weights = X_sim_frame$adm2_pop_weight,
    #                 level = .8)
    # bbin_LGM_simframe_res$betabinomial.spde.lgm.est.subdomain$method <-
    #   "Betabinomial SPDE LGM Sim. Frame"
    # holdout_i <- match(holdout_name,
    #                    bbin_LGM_simframe_res$betabinomial.spde.lgm.est.subdomain$domain)
    # res_list <- 
    #   c(res_list,
    #     list(bbin_LGM_simframe_res$betabinomial.spde.lgm.est.subdomain[holdout_i, ]))
    # 
    
    # nested version of admin2_mat
    admin2_mat_nested <- admin2_mat
    for (i in 1:nrow(admin2_mat)) {
      admin2_mat_nested[i, which(poly_adm2$NAME_1 != poly_adm2$NAME_1[i])] <- 0
      if (sum(admin2_mat_nested[i, ]) > 0) {
        admin2_mat_nested[i,] <- admin2_mat_nested[i,] / sum(admin2_mat_nested[i,])
      }
    }
    X_pop$admin1_name[X_pop$admin1_name == "North-Western"] <- "NorthWestern"
    
    holdout_des$variables$admin1_name[holdout_des$variables$admin1_name == "North-Western"] <- "NorthWestern"
    #### BINOMIAL NESTED ####
    print("NESTED MODELS")
    bin_nested_adm2_res <- 
      smoothUnitNested(
        hiv ~ admin1_name + urban + l1pntl*urban + l1paccess*urban + malaria*urban,
        domain = ~admin2_name, 
        design = holdout_des, 
        adj.mat = admin2_mat_nested,
        X.pop = X_pop,
        family = "binomial",
        X.pop.weights = X_pop$adm2_pop_weight,
        level = .8,
        return.samples = T
      )
    bin_nested_adm2_res$bym2.model.est$method <- "Binomial BYM2 Nested"
    holdout_i <- match(holdout_name,
                       bin_nested_adm2_res$bym2.model.est$domain)
    bym2_logit_samples <-
      SUMMER::logit(bin_nested_adm2_res$bym2.model.sample)[match(unique(adm2_est_table$domain), rownames(admin2_mat)), ]
    bin_nested_adm2_res$bym2.model.est$logit_mean <- 
      apply(bym2_logit_samples, 1, mean)
    bin_nested_adm2_res$bym2.model.est$logit_var <- 
      apply(bym2_logit_samples, 1, var)
    res_list <- 
      c(res_list,
        list(bin_nested_adm2_res$bym2.model.est[holdout_i, ]))
    
    
    bbin_nested_adm2_res <-
      fitBetabinomialBYM2(
        hiv ~ admin1_name + urban + l1pntl*urban + l1paccess*urban + malaria*urban,
        domain = ~admin2_name, 
        design = holdout_des, 
        adj.mat = admin2_mat_nested,
        X.pop = X_pop,
        X.pop.weights = X_pop$adm2_pop_weight,
        level = .8,
        return.samples = T
      )
    bbin_nested_adm2_res$betabinomial.bym2.model.est$method <- "Betabinomial BYM2 Nested"
    holdout_i <- match(holdout_name,
                       bbin_nested_adm2_res$betabinomial.bym2.model.est$domain)
    bym2_logit_samples <-
      SUMMER::logit(bbin_nested_adm2_res$betabinomial.bym2.model.sample)[match(unique(adm2_est_table$domain), rownames(admin2_mat)), ]
    bbin_nested_adm2_res$betabinomial.bym2.model.est$logit_mean <- 
      apply(bym2_logit_samples, 1, mean)
    bbin_nested_adm2_res$betabinomial.bym2.model.est$logit_var <- 
      apply(bym2_logit_samples, 1, var)
    res_list <- 
      c(res_list,
        list(bbin_nested_adm2_res$betabinomial.bym2.model.est[holdout_i, ]))
    
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

comp_table <- holdout_res |>
  group_by(method) |>
  summarize(mae = mean(abs_err),
            rmse = sqrt(mean(sq_err)),
            cov = mean(cov),
            log_cpo = mean(log_cpo),
            is = mean(is))

cov = mean(lower < direct.est & upper > direct.est),
log_cpo = mean(dnorm(SUMMER::logit(direct.est),
                     mean = SUMMER::logit(median), 
                     sd = sqrt(direct.est.var / direct.est^2 / (1-direct.est) ^2), log = T)),
is = mean((upper - lower) + 
            2 / 0.2 * ifelse(direct.est < lower, lower - direct.est, 0) +
            2 / 0.2 * ifelse(direct.est > upper, direct.est - upper, 0)))

saveRDS(holdout_res, "results/estimates/holdout_res.rds")

holdout_res |> 
  group_by(method) |> 
  summarize(rmse = sqrt(mean((median -direct.est) ^ 2)),
            mae = mean(abs(median - direct.est)),
            cov = mean(lower < direct.est & upper > direct.est))


holdout_res$method <- 
  rep(
    c(
      "Area level model: IID",
      "Area level model: BYM2",
      "Area level model: IID: No cov.",
      "Area level model: BYM2: No cov.",
      "Binomial BYM2",
      "Binomial BYM2: No cov.",
      "Lonobinomial BYM2",
      "Lonobinomial BYM2: No cov.",
      "Betabinomial BYM2",
      "Betabinomial BYM2: No cov.",
      "Binomial SPDE LGM",
      "Binomial SPDE LGM: No cov.",
      "Lonobinomial SPDE LGM",
      "Lonobinomial SPDE LGM: No cov.",
      "Betabinomial SPDE LGM",
      "Betabinomial SPDE LGM: No cov.",
      "Lonobinomial BYM2 Sim. Frame",
      "Betabinomial BYM2 Sim. Frame",
      "Lonobinomial SPDE LGM Sim. Frame",
      "Betabinomial SPDE LGM Sim. Frame"
    ), 
    22
  )



