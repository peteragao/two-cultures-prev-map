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


# LOAD RESULTS -----------------------------------------------------------------
adm1_alm_res <- readRDS("results/estimates/adm1_alm_res.rds")
adm2_alm_res <- readRDS("results/estimates/adm2_alm_res.rds")

adm1_cov_alm_res <- readRDS("results/estimates/adm1_cov_alm_res.rds")
adm2_cov_alm_res <- readRDS("results/estimates/adm2_cov_alm_res.rds")

bin_adm1_res <- readRDS("results/estimates/bin_adm1_res.rds")
bin_adm2_res <- readRDS("results/estimates/bin_adm2_res.rds")

bbin_adm1_res <- readRDS("results/estimates/bbin_adm1_res.rds")
bbin_adm2_res <- readRDS("results/estimates/bbin_adm2_res.rds")
bbin_no_cov_adm1_res <- readRDS("results/estimates/bbin_no_cov_adm1_res.rds")
bbin_no_cov_adm2_res <- readRDS("results/estimates/bbin_no_cov_adm2_res.rds")


bin_LGM_fit <- readRDS("results/estimates/bin_LGM_fit.rds")
bin_LGM_res <- readRDS("results/estimates/bin_LGM_res.rds")

bbin_LGM_fit <- readRDS("results/estimates/bbin_LGM_fit.rds")
bbin_LGM_res <- readRDS("results/estimates/bbin_LGM_res.rds")


bbin_LGM_no_cov_fit <- readRDS("results/estimates/bbin_LGM_no_cov_fit.rds")
bbin_LGM_no_cov_res <- readRDS("results/estimates/bbin_LGM_no_cov_res.rds")

adm1_est_table <- read.csv("results/estimates/adm1_est.csv")
adm2_est_table <- read.csv("results/estimates/adm2_est.csv")

selected_methods <- 
  c("Direct",
    "Fay-Herriot BYM2",
    "Betabinomial BYM2 No Cov.",
    "Betabinomial BYM2",
    "Betabinomial GRF No Cov.",
    "Betabinomial GRF")

adm1_est_table <- adm1_est_table |>
  filter(method %in% selected_methods) |>
  mutate(method = factor(method, levels = selected_methods)) 

adm2_est_table <- adm2_est_table |>
  filter(method %in% selected_methods) |>
  mutate(method = factor(method, levels = selected_methods)) 
# FIGURES ----------------------------------------------------------------------
my_theme <- 
  theme(legend.position="bottom",
        strip.text = element_text(size=20),
        legend.text = element_text(size=20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border=element_blank(),
        axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        legend.key.size = unit(1, "cm"),
        legend.title = element_text(size=20),
        legend.key.width = unit(1.5,"cm"))



#### TABLE 1: MODEL PARAMETERS ####

median_cl_size <- svy_dat |> 
  st_set_geometry(NULL) |>
  group_by(cluster) |>
  summarize(n = n()) |>
  summarize(n = median(n)) |>
  pull(n)

model_summaries <- data.frame(
  Method = 
    rep(c("Fay-Herriot BYM2 No Cov.",
          "Fay-Herriot BYM2",
          "Betabinomial BYM2 No Cov.",
          "Betabinomial BYM2",
          "Betabinomial GRF No Cov.",
          "Betabinomial GRF"), 2),
  `Interval Width` = 
    c(mean(adm1_alm_res$bym2.model.est$upper - 
             adm1_alm_res$bym2.model.est$lower),
      mean(adm1_cov_alm_res$bym2.model.est$upper - 
             adm1_cov_alm_res$bym2.model.est$lower),
      mean(bbin_no_cov_adm1_res$betabinomial.bym2.model.est$upper - 
             bbin_no_cov_adm1_res$betabinomial.bym2.model.est$lower),
      mean(bbin_adm1_res$betabinomial.bym2.model.est$upper - 
             bbin_adm1_res$betabinomial.bym2.model.est$lower),
      mean(bbin_LGM_no_cov_res$betabinomial.spde.lgm.est$upper - 
             bbin_LGM_no_cov_res$betabinomial.spde.lgm.est$lower),
      mean(bbin_LGM_res$betabinomial.spde.lgm.est$upper - 
             bbin_LGM_res$betabinomial.spde.lgm.est$lower),
      mean(adm2_alm_res$bym2.model.est$upper - 
             adm2_alm_res$bym2.model.est$lower),
      mean(adm2_cov_alm_res$bym2.model.est$upper - 
             adm2_cov_alm_res$bym2.model.est$lower),
      mean(bbin_no_cov_adm2_res$betabinomial.bym2.model.est$upper - 
             bbin_no_cov_adm2_res$betabinomial.bym2.model.est$lower),
      mean(bbin_adm2_res$betabinomial.bym2.model.est$upper - 
             bbin_adm2_res$betabinomial.bym2.model.est$lower),
      mean(bbin_LGM_no_cov_res$betabinomial.spde.lgm.est.subdomain$upper - 
             bbin_LGM_no_cov_res$betabinomial.spde.lgm.est.subdomain$lower),
      mean(bbin_LGM_res$betabinomial.spde.lgm.est.subdomain$upper - 
             bbin_LGM_res$betabinomial.spde.lgm.est.subdomain$lower)),
  Overdisp = 
    c(NA, NA,
      1 + (median_cl_size - 1) * 
        bbin_no_cov_adm1_res$betabinomial.bym2.model.fit$summary.hyperpar[1, 4],
      1 + (median_cl_size - 1) * 
        bbin_adm1_res$betabinomial.bym2.model.fit$summary.hyperpar[1, 4],
      1 + (median_cl_size - 1) * 
        bbin_LGM_no_cov_fit$summary.hyperpar[1, 4],
      1 + (median_cl_size - 1) * 
        bbin_LGM_fit$summary.hyperpar[1, 4],
      NA, NA,
      1 + (median_cl_size - 1) * 
        bbin_no_cov_adm2_res$betabinomial.bym2.model.fit$summary.hyperpar[1, 4],
      1 + (median_cl_size - 1) * 
        bbin_adm2_res$betabinomial.bym2.model.fit$summary.hyperpar[1, 4],
      1 + (median_cl_size - 1) * 
        bbin_LGM_no_cov_fit$summary.hyperpar[1, 4],
      1 + (median_cl_size - 1) * 
        bbin_LGM_fit$summary.hyperpar[1, 4]),
  `Spatial SD` = 
    c(1/sqrt(adm1_alm_res$bym2.model.fit$summary.hyperpar[1, 4]),
      1/sqrt(adm1_cov_alm_res$bym2.model.fit$summary.hyperpar[1, 4]),
      1/sqrt(bbin_no_cov_adm1_res$betabinomial.bym2.model.fit$summary.hyperpar[2, 4]),
      1/sqrt(bbin_adm1_res$betabinomial.bym2.model.fit$summary.hyperpar[2, 4]),
      bbin_LGM_no_cov_fit$summary.hyperpar[3, 4],
      bbin_LGM_fit$summary.hyperpar[3, 4],
      1/sqrt(adm2_alm_res$bym2.model.fit$summary.hyperpar[1, 4]),
      1/sqrt(adm2_cov_alm_res$bym2.model.fit$summary.hyperpar[1, 4]),
      1/sqrt(bbin_no_cov_adm2_res$betabinomial.bym2.model.fit$summary.hyperpar[2, 4]),
      1/sqrt(bbin_adm2_res$betabinomial.bym2.model.fit$summary.hyperpar[2, 4]),
      bbin_LGM_no_cov_fit$summary.hyperpar[3, 4],
      bbin_LGM_fit$summary.hyperpar[3, 4]
    )
)

# Update since the models have a mix of CI levels
getCI <- function(x, CI = 0.95){quantile(x, 1 - (1 - CI)/2) - quantile(x, (1 - CI)/2)}
# Admin1
model_summaries[, 2] <- c(
  mean(apply(adm1_alm_res$bym2.model.sample, 1, getCI)),
  mean(apply(adm1_cov_alm_res$bym2.model.sample, 1, getCI)),
  mean(apply(bbin_no_cov_adm1_res$betabinomial.bym2.model.sample, 1, getCI)),
  mean(apply(bbin_adm1_res$betabinomial.bym2.model.sample, 1, getCI)),
  mean(apply(bbin_LGM_no_cov_res$betabinomial.spde.lgm.sample, 1, getCI)),
  mean(apply(bbin_LGM_res$betabinomial.spde.lgm.sample, 1, getCI)),

  # Admin2,
  mean(apply(adm2_alm_res$bym2.model.sample, 1, getCI)),
  mean(apply(adm2_cov_alm_res$bym2.model.sample, 1, getCI)),
  mean(apply(bbin_no_cov_adm2_res$betabinomial.bym2.model.sample, 1, getCI)),
  mean(apply(bbin_adm2_res$betabinomial.bym2.model.sample, 1, getCI)),
  mean(apply(bbin_LGM_no_cov_res$betabinomial.spde.lgm.sample.subdomain, 1, getCI)),
  mean(apply(bbin_LGM_res$betabinomial.spde.lgm.sample.subdomain, 1, getCI))
)

model_summaries |>
  mutate(Interval.Width = Interval.Width * 100) |>
  knitr::kable( format = "latex", booktabs = T, linesep = "", digits = 2) 
  # |>
  # writeLines(paste0("results/figures/zmb_hiv_model_summaries.tex"))



#### ADMIN 1 ESTIMATES ####
plot_sf <- poly_adm1 |>
  left_join(adm1_est_table, by = "domain")
reordered_methods <- 
  c("Direct",
    "Betabinomial BYM2 No Cov.",
    "Betabinomial BYM2",
    "Fay-Herriot BYM2",
    "Betabinomial GRF No Cov.",
    "Betabinomial GRF")

plot_sf <- plot_sf |>
  mutate(method = factor(method, levels = reordered_methods)) |>
  mutate(CV = sqrt(var) / median * 100)
ggplot(plot_sf) + 
  geom_sf(aes(fill = median)) +
  facet_wrap(~method) +
  scale_fill_whitebox_c(palette = "viridi", direction = -1, limits = c(0, .33)) +
  theme_minimal()  + my_theme 
ggsave("results/figures/Zambia_adm1_est_map.png",
       width = 12.5, height = 8)
#### ADMIN 1 CVS ####
plot_sf <- plot_sf |>
  mutate(CV = sqrt(var) / median * 100)
ggplot(plot_sf) + 
  geom_sf(aes(fill = CV)) +
  facet_wrap(~method) +
  scale_fill_viridis_c(option = "magma", direction = -1, name = "CV (%)") +
  theme_minimal()  + my_theme 
ggsave("results/figures/Zambia_adm1_cv_map.png",
       width = 12.5, height = 8)

#### ADMIN 2 ESTIMATES ####
plot_sf <- poly_adm2 |>
  left_join(adm2_est_table, by = "domain")
reordered_methods <- 
  c("Direct",
    "Betabinomial BYM2 No Cov.",
    "Betabinomial BYM2",
    "Fay-Herriot BYM2",
    "Betabinomial GRF No Cov.",
    "Betabinomial GRF")

plot_sf <- plot_sf |>
  mutate(method = factor(method, levels = reordered_methods)) |>
  mutate(CV = sqrt(var) / median * 100)
ggplot(plot_sf) + 
  geom_sf(aes(fill = median)) +
  facet_wrap(~method,nrow = 2) +
  scale_fill_whitebox_c(palette = "viridi", direction = -1, limits = c(0, .33)) +
  theme_minimal()  + my_theme
ggsave("results/figures/Zambia_adm2_est_map.png",
       width = 12.5, height = 8)
#### ADMIN 2 CVs ####
plot_sf <- plot_sf |>
  mutate(CV = sqrt(var) / median * 100)
ggplot(plot_sf) + 
  geom_sf(aes(fill = CV)) +
  facet_wrap(~method, nrow = 2) +
  scale_fill_viridis_c(option = "magma", direction = -1, name = "CV (%)", na.value="white") +
  theme_minimal()  + my_theme 
ggsave("results/figures/Zambia_adm2_cv_map.png",
       width = 12.5, height = 8)
#### TABLE 2: CV RESULTS ####
adm1_holdout_res <- readRDS("results/estimates/adm1_holdout_res.rds")

selected_methods <- 
  c("Direct",
    "Fay-Herriot BYM2",
    "Betabinomial BYM2: No cov.",
    "Betabinomial BYM2",
    "Betabinomial GRF: No cov.",
    "Betabinomial GRF")
adm1_holdout_res <- adm1_holdout_res |>
  filter(method %in% selected_methods) |>
  mutate(method = factor(method, levels = selected_methods)) |>
  left_join(adm1_alm_res$direct.est |>
              select(domain, median, var) |>
              rename(direct.est = median, 
                     direct.est.var = var),
            by = "domain") |>
  filter(!is.na(direct.est)) |>
  mutate(logit_direct_est = SUMMER::logit(direct.est),
         logit_direct_est_var = direct.est.var / direct.est^2 / (1-direct.est) ^2) 



adm1_holdout_res <- adm1_holdout_res |>
  mutate(logit_upper = logit_mean + qnorm(.9) * sqrt(logit_direct_est_var + logit_var),
         logit_lower = logit_mean - qnorm(.9) * sqrt(logit_direct_est_var + logit_var),
         sq_err = (logit_mean - logit_direct_est)^2,
         abs_err = abs(logit_mean - logit_direct_est),
         cov = logit_lower < logit_direct_est & logit_upper > logit_direct_est,
         log_cpo = dnorm(logit_direct_est, mean = logit_mean, sd = sqrt(logit_direct_est_var + logit_var), log = T),
         is = logit_upper - logit_lower + 
           2 / 0.2 * ifelse(logit_direct_est < logit_lower, logit_lower - logit_direct_est, 0) +
           2 / 0.2 * ifelse(logit_direct_est > logit_upper, logit_direct_est - logit_upper, 0))
comp_table <- adm1_holdout_res |>
  group_by(method) |>
  summarize(mae = mean(abs_err),
            rmse = sqrt(mean(sq_err)),
            cov = mean(cov),
            log_cpo = mean(log_cpo),
            is = mean(is),
            int_length = mean(logit_upper - logit_lower))
present_comp_table <- comp_table |>
  mutate(mae = round(mae * 100, 2),
         rmse = round(rmse * 100, 2),
         cov = round(cov * 100),
         log_cpo = round(log_cpo - max(log_cpo), 2),
         is = round(is, 2),
         int_length = round(int_length * 100, 2))
full_res_list <- list()
for (holdout_name in poly_adm2$NAME_2) {
  if (file.exists(paste0("results/estimates/adm2_holdout/holdout_res_", holdout_name, ".rds"))) {
    print(holdout_name)
    holdout_res <- readRDS(paste0("results/estimates/adm2_holdout/holdout_res_", holdout_name, ".rds"))
    full_res_list <- c(full_res_list, list(holdout_res))
  }
}
adm2_holdout_res <- bind_rows(full_res_list)
adm2_holdout_res <- adm2_holdout_res |>
  filter(method %in% selected_methods) |>
  mutate(method = factor(method, levels = selected_methods)) |>
  left_join(adm2_alm_res$direct.est |>
              select(domain, median, var) |>
              rename(direct.est = median, 
                     direct.est.var = var),
            by = "domain") |>
  filter(!is.na(direct.est)) |>
  mutate(logit_direct_est = SUMMER::logit(direct.est),
         logit_direct_est_var = direct.est.var / direct.est^2 / (1-direct.est) ^2) 
adm2_holdout_res <- adm2_holdout_res |>
  mutate(logit_upper = logit_mean + qnorm(.9) * sqrt(logit_direct_est_var + logit_var),
         logit_lower = logit_mean - qnorm(.9) * sqrt(logit_direct_est_var + logit_var),
         sq_err = (logit_mean - logit_direct_est)^2,
         abs_err = abs(logit_mean - logit_direct_est),
         cov = logit_lower < logit_direct_est & logit_upper > logit_direct_est,
         log_cpo = dnorm(logit_direct_est, mean = logit_mean, sd = sqrt(logit_direct_est_var + logit_var), log = T),
         is = logit_upper - logit_lower + 
           2 / 0.2 * ifelse(logit_direct_est < logit_lower, logit_lower - logit_direct_est, 0) +
           2 / 0.2 * ifelse(logit_direct_est > logit_upper, logit_direct_est - logit_upper, 0))
comp_table <- adm2_holdout_res |>
  group_by(method) |>
  summarize(mae = mean(abs_err),
            rmse = sqrt(mean(sq_err)),
            cov = mean(cov),
            log_cpo = mean(log_cpo),
            is = mean(is),
            int_length = mean(logit_upper - logit_lower))
present_comp_table <- comp_table |>
  mutate(mae = round(mae * 100, 2),
         rmse = round(rmse * 100, 2),
         cov = round(cov * 100),
         log_cpo = round(log_cpo - max(log_cpo), 2),
         is = round(is, 2),
         int_length = round(int_length * 100, 2))

ggplot(adm2_holdout_res, aes(x = logit_direct_est, y = logit_mean, color = method)) + 
  geom_point()+
  facet_wrap(~method)
ggplot(adm1_holdout_res, aes(x = logit_direct_est, y = logit_mean, color = method)) + 
  geom_point()+
  facet_wrap(~method) +
  geom_abline(slope =1 )


#### ******** APPENDIX ******** ####
#### FIG S1 COVARIATE MAPS ####
cen_pop <- rast("data/Zambia/Population/zmb_ppp_2010.tif")
cen_pop_1km <- aggregate(cen_pop, fact = 10, fun = "sum", na.rm = T)
covs <- c("urban", "access", "l1pntl", "malaria")
X_pop_sf <- st_as_sf(X_pop, coords = c("LONGNUM", "LATNUM"), crs = "EPSG:4326")

rast_list <- lapply(covs, function(x) {
  rasterize(X_pop_sf, cen_pop_1km,
            field = x)
})
cov_rast <- do.call("c", rast_list)
names(cov_rast) <- covs
urban_gg <- ggplot() +
  geom_spatraster(data = cov_rast, aes(fill = urban)) +
  scale_fill_whitebox_c(palette = "viridi") +
  labs(fill = "", title = "Urban") +
  theme_minimal(base_size = 22)  + my_theme
access_gg <- ggplot() +
  geom_spatraster(data = cov_rast, aes(fill = access)) +
  scale_fill_whitebox_c(palette = "viridi") +
  labs(fill = "", title = "Access") +
  theme_minimal(base_size = 22)  + my_theme
ntl_gg <- ggplot() +
  geom_spatraster(data = cov_rast, aes(fill = l1pntl)) +
  scale_fill_whitebox_c(palette = "viridi") +
  labs(fill = "",
       title = "Night time lights") +
  theme_minimal(base_size = 22)  + my_theme
malaria_gg <- ggplot() +
  geom_spatraster(data = cov_rast, aes(fill = malaria)) +
  scale_fill_whitebox_c(palette = "viridi") +
  labs(fill = "", title = "Malaria prevalence") +
  theme_minimal(base_size = 22)  + my_theme
library(patchwork)
cov_gg <- (urban_gg + access_gg) / (ntl_gg + malaria_gg)
ggsave("results/figures/Zambia_covariate_map.png",
       cov_gg,
       width = 10, height = 11.4)
#### TABLE S2 CROSS-VALIDATION ####
#### FIG S3 INTERVALS SCATTER ####

ggplot(adm2_est_table |>
         filter(!(method == "Fay-Herriot IID")),
       aes(x = direct.est, y = mean, color = method)) +
  geom_point(size = 3, show.legend = F) +
  geom_abline(slope = 1) +
  geom_linerange(aes(ymin = lower, ymax = upper), show.legend = F) + 
  facet_wrap(~method, nrow =3) + 
  xlab("Direct Estimate") +
  ylab("Model Estimate") +
  scale_color_manual(values = RColorBrewer::brewer.pal(6, "Set1")) +
  theme_minimal(base_size = 24) + 
  theme(
    legend.position="bottom",
    panel.grid.minor = element_blank(),
    legend.key.size = unit(.5, "cm"),
    legend.key.width = unit(.5,"cm")
  )
ggsave("results/figures/Zambia_adm2_scatter_bars.png",
       width = 8, height = 12)


#### ADMIN 1 NATIONAL AGGREGATES ####
adm1_to_natl_prop <- X_pop |> 
  mutate(natl_f1549_pop = sum(Population)) |>
  group_by(admin1_name) |> 
  mutate(admin1_f1549_pop_prop = sum(Population) / natl_f1549_pop) |>
  ungroup() |>
  dplyr::select(admin1_name, admin1_f1549_pop_prop) |>
  unique() |> 
  arrange(match(admin1_name, bin_adm1_res$bym2.model.est$domain))
adm1_alm_iid_sample <-
  matrix(
    rnorm(10*250, 
          mean = adm1_alm_res$iid.model.est$mean, 
          sd = sqrt(adm1_alm_res$iid.model.est$var)),
    ncol = 250
  )
adm1_alm_bym2_sample <-
  matrix(
    rnorm(10*250, 
          mean = adm1_alm_res$bym2.model.est$mean, 
          sd = sqrt(adm1_alm_res$bym2.model.est$var)),
    ncol = 250
  )
adm1_cov_alm_bym2_sample <-
  matrix(
    rnorm(10*250, 
          mean = adm1_cov_alm_res$bym2.model.est$mean, 
          sd = sqrt(adm1_cov_alm_res$bym2.model.est$var)),
    ncol = 250
  )
adm1_alm_iid_sample <- t(adm1_alm_iid_sample) %*% adm1_to_natl_prop$admin1_f1549_pop_prop
adm1_alm_bym2_sample <- t(adm1_alm_bym2_sample) %*% adm1_to_natl_prop$admin1_f1549_pop_prop
adm1_cov_alm_bym2_sample <- t(adm1_cov_alm_bym2_sample) %*% adm1_to_natl_prop$admin1_f1549_pop_prop
bin_adm1_natl_sample <- 
  t(bin_adm1_res$bym2.model.sample)[, match(adm1_to_natl_prop$admin1_name, rownames(admin1_mat))] %*%
  adm1_to_natl_prop$admin1_f1549_pop_prop 
bbin_no_cov_adm1_natl_sample <- 
  t(bbin_no_cov_adm1_res$betabinomial.bym2.model.sample)[, match(adm1_to_natl_prop$admin1_name, rownames(admin1_mat))] %*%
  adm1_to_natl_prop$admin1_f1549_pop_prop 
bbin_adm1_natl_sample <- 
  t(bbin_adm1_res$betabinomial.bym2.model.sample)[, match(adm1_to_natl_prop$admin1_name, rownames(admin1_mat))] %*%
  adm1_to_natl_prop$admin1_f1549_pop_prop 
bin_LGM_natl_sample <- 
  t(bin_LGM_res$binomial.spde.lgm.sample)[, match(adm1_to_natl_prop$admin1_name, rownames(admin1_mat))] %*%
  adm1_to_natl_prop$admin1_f1549_pop_prop 
bbin_LGM_no_cov_natl_sample <- 
  t(bbin_LGM_no_cov_res$betabinomial.spde.lgm.sample)[, match(adm1_to_natl_prop$admin1_name, rownames(admin1_mat))] %*%
  adm1_to_natl_prop$admin1_f1549_pop_prop 
bbin_LGM_natl_sample <- 
  t(bbin_LGM_res$betabinomial.spde.lgm.sample)[, match(adm1_to_natl_prop$admin1_name, rownames(admin1_mat))] %*%
  adm1_to_natl_prop$admin1_f1549_pop_prop 

natl_comp <- data.frame(
  method =
    c("Direct (National)",
      # "Fay-Herriot IID",
      "Fay-Herriot BYM2: No Cov.",
      "Fay-Herriot BYM2",
      # "Binomial BYM2",
      "Betabinomial BYM2: No Cov.",
      "Betabinomial BYM2",
      # "Binomial GRF",
      "Betabinomial GRF: No Cov.",
      "Betabinomial GRF"),
  est = 
    c(
      svymean(~hiv, sample_des),
      # median(adm1_alm_iid_sample),
      median(adm1_alm_bym2_sample),
      median(adm1_cov_alm_bym2_sample),
      # median(bin_adm1_natl_sample),
      median(bbin_no_cov_adm1_natl_sample),
      median(bbin_adm1_natl_sample),
      # median(bin_LGM_natl_sample),
      median(bbin_LGM_no_cov_natl_sample),
      median(bbin_LGM_natl_sample)
    ),
  se = 
    c(
      sqrt(vcov(svymean(~hiv, sample_des))),
      # sd(adm1_alm_iid_sample),
      sd(adm1_alm_bym2_sample),
      sd(adm1_cov_alm_bym2_sample),
      # sd(bin_adm1_natl_sample),
      sd(bbin_no_cov_adm1_natl_sample),
      sd(bbin_adm1_natl_sample),
      # sd(bin_LGM_natl_sample),
      sd(bbin_LGM_no_cov_natl_sample),
      sd(bbin_LGM_natl_sample)
    )
) |>
  mutate(lower = est - qnorm(.975) * se,
         upper = est + qnorm(.975) * se) |>
  mutate(method =factor(method,
                        levels =   c("Direct (National)",
                                     # "Fay-Herriot IID",
                                     "Fay-Herriot BYM2: No Cov.",
                                     "Fay-Herriot BYM2",
                                     # "Binomial BYM2",
                                     "Betabinomial BYM2: No Cov.",
                                     "Betabinomial BYM2",
                                     # "Binomial GRF",
                                     "Betabinomial GRF: No Cov.",
                                     "Betabinomial GRF")))
g <- ggplot(natl_comp, aes(x = method, y = est, color = method == "Direct (National)")) +
  geom_linerange(aes(x = method, ymin = lower, ymax = upper)) +
  geom_point() +
  xlab("National estimate") +
  ylab("Method") +
  scale_color_manual(values = c("black", "red")) + 
  ggtitle("Aggregated Admin-1 level model-based estimates\nvs. national direct estimate") +
  theme_minimal(base_size = 20) +
  ggplot2::theme(legend.position = "none",
                 axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust=1))
ggsave(g, file = paste0("results/figures/Zambia_adm1_aggregated_vs_natl_direct.png"),
       width = 9, height = 7)

#### ADMIN1 ESTIMATES TABLE ####
out_table <- adm1_est_table |>
  filter(method %in% selected_methods) |>
  filter(complete.cases(median)) |>
  mutate(est = paste0(round(median, 2), " (", round(lower, 2),", ", round(upper, 2), ")")) |>
  dplyr::select(est, median, method, domain) |>
  pivot_wider(values_from = c("median", "est"),
              names_from = method,
              names_glue = "{method}_{.value}") |>
  arrange(Direct_median) |>
  dplyr::select(domain, contains("est")) 
names(out_table) <- 
  c("Province",
    stringr::str_replace(names(out_table)[-1], "_est", ""))

out_table |>
  knitr::kable( format = "latex", booktabs = T, linesep = "", digits = 3) |>
  writeLines(paste0("results/figures/zmb_hiv_adm1_est_table.tex"))




cv_table <- adm2_est_table |>
  filter(method == "Direct") |>
  mutate(CV = sqrt(var) / median * 100) |>
  dplyr::select(domain, median, CV) |>
  arrange(CV) 

cv_table |>
  knitr::kable( format = "latex", booktabs = T, linesep = "", digits = 3) |>
  writeLines(paste0("results/figures/zmb_hiv_adm2_cv_table.tex"))

#### ADMIN 2 NATIONAL AGGREGATES ####
adm2_to_natl_prop <- X_pop |> 
  mutate(natl_f1549_pop = sum(Population)) |>
  group_by(admin2_name) |> 
  mutate(admin2_f1549_pop_prop = sum(Population) / natl_f1549_pop) |>
  ungroup() |>
  dplyr::select(admin2_name, admin2_f1549_pop_prop) |>
  unique() |> 
  arrange(match(admin2_name, bin_adm2_res$bym2.model.est$domain))
adm2_alm_iid_sample <-
  matrix(
    rnorm(115*250, 
          mean = adm2_alm_res$iid.model.est$mean, 
          sd = sqrt(adm2_alm_res$iid.model.est$var)),
    ncol = 250
  )
adm2_alm_bym2_sample <-
  matrix(
    rnorm(115*250, 
          mean = adm2_alm_res$bym2.model.est$mean, 
          sd = sqrt(adm2_alm_res$bym2.model.est$var)),
    ncol = 250
  )
adm2_cov_alm_bym2_sample <-
  matrix(
    rnorm(115*250, 
          mean = adm2_cov_alm_res$bym2.model.est$mean, 
          sd = sqrt(adm2_cov_alm_res$bym2.model.est$var)),
    ncol = 250
  )
adm2_alm_iid_sample <- t(adm2_alm_iid_sample) %*% adm2_to_natl_prop$admin2_f1549_pop_prop
adm2_alm_bym2_sample <- t(adm2_alm_bym2_sample) %*% adm2_to_natl_prop$admin2_f1549_pop_prop
adm2_cov_alm_bym2_sample <- t(adm2_cov_alm_bym2_sample) %*% adm2_to_natl_prop$admin2_f1549_pop_prop
bin_adm2_natl_sample <- 
  t(bin_adm2_res$bym2.model.sample)[, match(adm2_to_natl_prop$admin2_name, rownames(admin2_mat))] %*%
  adm2_to_natl_prop$admin2_f1549_pop_prop 
bbin_no_cov_adm2_natl_sample <- 
  t(bbin_no_cov_adm2_res$betabinomial.bym2.model.sample)[, match(adm2_to_natl_prop$admin2_name, rownames(admin2_mat))] %*%
  adm2_to_natl_prop$admin2_f1549_pop_prop 
bbin_adm2_natl_sample <- 
  t(bbin_adm2_res$betabinomial.bym2.model.sample)[, match(adm2_to_natl_prop$admin2_name, rownames(admin2_mat))] %*%
  adm2_to_natl_prop$admin2_f1549_pop_prop 
bin_LGM_natl_sample <- 
  t(bin_LGM_res$binomial.spde.lgm.sample.subdomain)[, match(adm2_to_natl_prop$admin2_name, rownames(admin2_mat))] %*%
  adm2_to_natl_prop$admin2_f1549_pop_prop 
bbin_LGM_natl_sample <- 
  t(bbin_LGM_res$betabinomial.spde.lgm.sample.subdomain)[, match(adm2_to_natl_prop$admin2_name, rownames(admin2_mat))] %*%
  adm2_to_natl_prop$admin2_f1549_pop_prop 

natl_comp <- data.frame(
  method =
    c("Direct (National)",
      # "Fay-Herriot IID",
      "Fay-Herriot BYM2: No Cov.",
      "Fay-Herriot BYM2",
      # "Binomial BYM2",
      "Betabinomial BYM2: No Cov.",
      "Betabinomial BYM2",
      # "Binomial GRF",
      "Betabinomial GRF: No Cov.",
      "Betabinomial GRF"),
  est = 
    c(
      svymean(~hiv, sample_des),
      # median(adm1_alm_iid_sample),
      median(adm2_alm_bym2_sample),
      median(adm2_cov_alm_bym2_sample),
      # median(bin_adm1_natl_sample),
      median(bbin_no_cov_adm2_natl_sample),
      median(bbin_adm2_natl_sample),
      # median(bin_LGM_natl_sample),
      median(bbin_LGM_no_cov_natl_sample),
      median(bbin_LGM_natl_sample)
    ),
  se = 
    c(
      sqrt(vcov(svymean(~hiv, sample_des))),
      # sd(adm1_alm_iid_sample),
      sd(adm2_alm_bym2_sample),
      sd(adm2_cov_alm_bym2_sample),
      # sd(bin_adm1_natl_sample),
      sd(bbin_no_cov_adm2_natl_sample),
      sd(bbin_adm2_natl_sample),
      # sd(bin_LGM_natl_sample),
      sd(bbin_LGM_no_cov_natl_sample),
      sd(bbin_LGM_natl_sample)
    )
) |>
  mutate(lower = est - qnorm(.975) * se,
         upper = est + qnorm(.975) * se) |>
    mutate(method =factor(method,
                        levels =   c("Direct (National)",
                                     # "Fay-Herriot IID",
                                     "Fay-Herriot BYM2: No Cov.",
                                     "Fay-Herriot BYM2",
                                     # "Binomial BYM2",
                                     "Betabinomial BYM2: No Cov.",
                                     "Betabinomial BYM2",
                                     # "Binomial GRF",
                                     "Betabinomial GRF: No Cov.",
                                     "Betabinomial GRF")))
ggplot(natl_comp, aes(x = method, y = est, color = method == "Direct (National)")) +
  geom_linerange(aes(x = method, ymin = lower, ymax = upper)) +
  geom_point() +
  xlab("National estimate") +
  ylab("Method") +
  scale_color_manual(values = c("black", "red")) + 
  ggtitle("Aggregated Admin-2 level model-based estimates\nvs. national direct estimate") +
  theme_minimal(base_size = 20) +
  ggplot2::theme(legend.position = "none",
                 axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust=1))
ggsave(paste0("results/figures/Zambia_adm2_aggregated_vs_natl_direct.png"),
       width = 9, height = 7)

#### ADMIN 2 ESTIMATES TABLE ####

out_table <- adm2_est_table |>
  filter(method %in% selected_methods) |>
  filter(complete.cases(median)) |>
  mutate(est = paste0(round(median, 2), " (", round(lower, 2),", ", round(upper, 2), ")")) |>
  dplyr::select(est, median, method, domain) |>
  pivot_wider(values_from = c("median", "est"),
              names_from = method,
              names_glue = "{method}_{.value}") |>
  arrange(Direct_median) |>
  dplyr::select(domain, contains("est")) 
names(out_table) <- 
  c("District",
    stringr::str_replace(names(out_table)[-1], "_est", ""))
out_table |>
  knitr::kable( format = "latex", booktabs = T, linesep = "", digits = 3) |>
  writeLines(paste0("results/figures/zmb_hiv_adm2_est_table.tex"))

#### ******** ADDITIONAL FIGURES ******** ####
#### ADDITIONAL MODEL PARAMETERS ####
adm1_alm_res <- readRDS("results/estimates/adm1_alm_res.rds")
adm2_alm_res <- readRDS("results/estimates/adm2_alm_res.rds")

bin_adm1_res <- readRDS("results/estimates/bin_adm1_res.rds")
bin_adm2_res <- readRDS("results/estimates/bin_adm2_res.rds")

bbin_adm1_res <- readRDS("results/estimates/bbin_adm1_res.rds")
bbin_adm2_res <- readRDS("results/estimates/bbin_adm2_res.rds")

bin_LGM_fit <- readRDS("results/estimates/bin_LGM_fit.rds")
bbin_LGM_fit <- readRDS("results/estimates/bbin_LGM_fit.rds")

unit_level_cov_adm1 <- bind_rows(
  list("Betabinomial BYM2 No Cov." =  bbin_no_cov_adm1_res$betabinomial.bym2.model.fit$summary.fixed,
       "Betabinomial BYM2" =  bbin_adm1_res$betabinomial.bym2.model.fit$summary.fixed,
       "Betabinomial GRF No Cov." = bbin_LGM_no_cov_fit$summary.fixed,
       "Betabinomial GRF" = bbin_LGM_fit$summary.fixed) |>
    map(function(x) tibble::rownames_to_column(x, "Variable")),
  .id = "method"
) |>
  mutate(Estimate = paste0(round(mean, 2), " (", round(sd, 2), ")")) |>
  mutate(Variable = stringr::str_replace(Variable, "intercept", "(Intercept)")) |>
  select(method, Variable, Estimate) |>
  pivot_wider(names_from = "method", values_from = "Estimate") |>
  knitr::kable( format = "latex", booktabs = T, linesep = "", digits = 3,
                caption = "Fixed effects estimates for unit-level covariate models (Admin-1)") |>
  writeLines(paste0("results/figures/zmb_hiv_adm1_unit_cov.tex"))

area_level_cov_adm1 <- bind_rows(
  list("Fay-Herriot BYM2" = adm1_alm_res$bym2.model.fit$summary.fixed) |>
    map(function(x) tibble::rownames_to_column(x, "Variable")),
  .id = "method"
) |>
  mutate(Estimate = paste0(round(mean, 2), " (", round(sd, 2), ")")) |>
  mutate(Variable = stringr::str_replace(Variable, "intercept", "(Intercept)")) |>
  select(method, Variable, Estimate) |>
  pivot_wider(names_from = "method", values_from = "Estimate") |>
  knitr::kable( format = "latex", booktabs = T, linesep = "", digits = 3,
                caption = "Fixed effects estimates for area-level covariate models (Admin-1)") |>
  writeLines(paste0("results/figures/zmb_hiv_adm1_area_cov.tex"))

unit_level_cov_adm2 <- bind_rows(
  list("Betabinomial BYM2 No Cov." =  bbin_no_cov_adm2_res$betabinomial.bym2.model.fit$summary.fixed,
       "Betabinomial BYM2" =  bbin_adm2_res$betabinomial.bym2.model.fit$summary.fixed,
       "Betabinomial GRF No Cov." = bbin_LGM_no_cov_fit$summary.fixed,
       "Betabinomial GRF" = bbin_LGM_fit$summary.fixed) |>
    map(function(x) tibble::rownames_to_column(x, "Variable")),
  .id = "method"
) |>
  mutate(Estimate = paste0(round(mean, 2), " (", round(sd, 2), ")")) |>
  mutate(Variable = stringr::str_replace(Variable, "intercept", "(Intercept)")) |>
  select(method, Variable, Estimate) |>
  pivot_wider(names_from = "method", values_from = "Estimate") |>
  knitr::kable(format = "latex", booktabs = T, linesep = "", digits = 3,
               caption = "Fixed effects estimates for unit-level covariate models (Admin-2)") |>
  writeLines(paste0("results/figures/zmb_hiv_adm2_unit_cov.tex"))

area_level_cov_adm2 <- bind_rows(
  list("Fay-Herriot BYM2" = adm2_alm_res$bym2.model.fit$summary.fixed) |>
    map(function(x) tibble::rownames_to_column(x, "Variable")),
  .id = "method"
) |>
  mutate(Estimate = paste0(round(mean, 2), " (", round(sd, 2), ")")) |>
  mutate(Variable = stringr::str_replace(Variable, "intercept", "(Intercept)")) |>
  select(method, Variable, Estimate) |>
  pivot_wider(names_from = "method", values_from = "Estimate") |>
  knitr::kable( format = "latex", booktabs = T, linesep = "", digits = 3,
                caption = "Fixed effects estimates for area-level covariate models (Admin-2)") |>
  writeLines(paste0("results/figures/zmb_hiv_adm2_area_cov.tex"))
hyperpar_list <- list(
  "Fay-Herriot BYM2 Admin-1" = adm1_alm_res$bym2.model.fit$summary.hyperpar,
  "Fay-Herriot BYM2 Admin-2" = adm2_alm_res$bym2.model.fit$summary.hyperpar,
  
  "Betabinomial BYM2 Admin-1 No Cov." = bbin_no_cov_adm1_res$betabinomial.bym2.model.fit$summary.hyperpar,
  "Betabinomial BYM2 Admin-1" = bbin_adm1_res$betabinomial.bym2.model.fit$summary.hyperpar,
  
  "Betabinomial BYM2 Admin-2 No Cov." = bbin_adm2_res$betabinomial.bym2.model.fit$summary.hyperpar,
  "Betabinomial BYM2 Admin-2" = bbin_no_cov_adm2_res$betabinomial.bym2.model.fit$summary.hyperpar,
  
  "Betabinomial GRF No Cov." = bbin_LGM_no_cov_fit$summary.hyperpar,
  "Betabinomial GRF" = bbin_LGM_fit$summary.hyperpar
)


f <- file("results/figures/zmb_hiv_hyperpar.tex", open = "wt")
for (i in 1:length(hyperpar_list)) {
  hyperpar_list[[i]] |> tibble::rownames_to_column("Variable") |>
    mutate(Estimate = paste0(round(mean, 2), " (", round(sd, 2), ")")) |>
    select(Variable, Estimate) |>
    knitr::kable(format = "latex", booktabs = T, linesep = "", digits = 3,
                 caption = paste0("Hyperparameter estimates for ", names(hyperpar_list)[i], " model")) |>
    write(f, append=TRUE)
}
close(f)

#### ADMIN 1 VS ADMIN 2 VARIABILITY ####
combined_table <- bind_rows(adm1_est_table |> mutate(level = "Admin-1"),
                            adm2_est_table |> mutate(level = "Admin-2"))
ggplot(combined_table, aes(x = median, y = level, color = method))+
  facet_wrap(~method) + 
  geom_boxplot(show.legend = F) +
  theme_minimal()  +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.position="bottom",
    panel.grid.minor = element_blank(),
    legend.key.size = unit(.5, "cm"),
    legend.text = element_text(size=14),
    legend.key.width = unit(.5,"cm"))
ggsave("results/figures/Zambia_adm1_vs_adm2_boxplots.png",
       width = 10, height = 8)

ggplot(combined_table, aes(x = median, y = level, color = method))+
  facet_wrap(~method) + 
  geom_violin(show.legend = F) +
  theme_minimal()  +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.position="bottom",
    panel.grid.minor = element_blank(),
    legend.key.size = unit(.5, "cm"),
    legend.text = element_text(size=14),
    legend.key.width = unit(.5,"cm"))
ggsave("results/figures/Zambia_adm1_vs_adm2_violins.png",
       width = 10, height = 8)


