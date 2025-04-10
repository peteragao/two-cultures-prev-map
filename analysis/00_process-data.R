
# 00_process-data.R
# This file loads and cleans covariate data and geographic boundaries.            
library(dplyr)
library(tidyr)
library(terra)
library(sf)
library(ggplot2)
library(tidyterra)
library(patchwork)
library(spdep)
country <- "Zambia"
survey_year <- 2018 
source("analysis/functions.R")
#### GADM POLYGON DATA ####
# GADM version 4.1 is used 
gadm_abbrev <- "ZMB"
poly_path <- "data/Zambia/GADM/gadm41_ZMB_shp"
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

# In 2018 Zambia redistricting converted the 72 districts into 116 districts
# Below, we visualize the new and old boundaries and create a table to convert
# the new districts to the old districts.
poly_path_old <- "data/Zambia/GADM/gadm36_ZMB_shp"
poly_layer_adm2_old <- paste('gadm36', gadm_abbrev,
                             '2', sep = "_")
poly_adm2_old <-  st_read(dsn = poly_path_old,
                          layer = as.character(poly_layer_adm2_old)) 

ggplot(poly_adm2) + 
  geom_sf(data = poly_adm2_old, color = "red", alpha = .5) + 
  geom_sf(data = poly_adm2, color = "blue", alpha = .5) + 
  theme_minimal()

district_table <- poly_adm2_old |>
  select(NAME_2) |>
  rename(District_2010 = NAME_2) |>
  st_intersection(poly_adm2 |>
                    select(NAME_2) |>
                    rename(District_2018 = NAME_2)) 
district_table <- district_table |>
  mutate(area = st_area(district_table)) |>
  group_by(District_2018) |>
  filter(area == max(area)) |>
  ungroup()




# POPULATION ESTIMATES ---------------------------------------------------------
load_pop <- function() {
  cen_pop <- rast("data/Zambia/Population/zmb_ppp_2010.tif")
  yrs <- seq(15, 49, by = 5)
  svy_f1549_pop_list <- 
    lapply(yrs, 
           function(yr) rast(paste0("data/Zambia/Population/zmb_f_",
                                    yr, "_2018.tif"))
    )
  svy_f1549_pop <- rast(svy_f1549_pop_list)
  svy_f1549_pop <- sum(svy_f1549_pop)
  
  cen_pop_1km <- aggregate(cen_pop, fact = 10, fun = "sum", na.rm = T)
  svy_f1549_pop_1km <- aggregate(svy_f1549_pop, fact = 10, fun = "sum", na.rm = T)
  
  pop_stk <- c(cen_pop_1km, svy_f1549_pop_1km)
  names(pop_stk) <- c("cen_pop", "svy_f1549_pop")
  adm1_pix_list <- list()
  adm1_pop_totals_list <- list()
  # make pixel table
  pop_pix <- st_join(st_as_sf(as.points(pop_stk)),
                     poly_adm2[, c("NAME_1", "NAME_2")],
                     join = st_nearest_feature)
  pop_pix <- rename(pop_pix, "admin1_name" = "NAME_1")
  pop_pix <- rename(pop_pix, "admin2_name" = "NAME_2")
  pop_pix$urban <- NA
  urban_rast <- rast("data/Zambia/Covariates/Zambia_UR_ind_surf.tif")
  pop_pix$urban <- extract(urban_rast, pop_pix)[,2]
  pop_pix <- pop_pix |> mutate(urban = tidyr::replace_na(urban, 0))
  adm2_pop_totals <- pop_pix |>
    st_set_geometry(NULL) |>
    group_by(admin1_name, admin2_name) |>
    summarize(cen_pop = sum(cen_pop),
              svy_f1549_pop = sum(svy_f1549_pop)) |>
    ungroup()
  adm2_urban_pop_totals <- pop_pix |>
    st_set_geometry(NULL) |>
    group_by(admin1_name, admin2_name, urban) |>
    summarize(cen_pop = sum(cen_pop),
              svy_f1549_pop = sum(svy_f1549_pop)) |>
    ungroup()
  adm2_urban_pop_totals <- subset(adm2_urban_pop_totals, urban == T)
  adm2_urban_pop_totals <- dplyr::select(adm2_urban_pop_totals, -urban)
  colnames(adm2_urban_pop_totals)[3:4] <- 
    paste0(colnames(adm2_urban_pop_totals)[3:4], "_urb")
  adm2_pop_totals <- merge(adm2_pop_totals, 
                           adm2_urban_pop_totals[, ], all.x = T)
  adm2_pop_totals[is.na(adm2_pop_totals)] <- 0
  saveRDS(pop_pix, 
          paste0("data/Zambia/zmb_pix_tbl.rds"))
  saveRDS(adm2_pop_totals, 
          paste0("data/Zambia/zmb_pop_totals.rds"))
}
pop_pix <- load_pop()
pop_pix <- readRDS(paste0("data/Zambia/zmb_pix_tbl.rds"))

# GEOSPATIAL COVARIATES --------------------------------------------------------


# Access (Travel time to cities in 2015 (Weiss et al. 2015))
pop_pix$access <-
  get_cov("data/Zambia/Covariates/2015_accessibility_to_cities_v1.0.tif", pop_pix)

# Night time lights (VIIRS (Worldpop 2016))
pop_pix$ntl <-
  get_cov("data/Zambia/Covariates/zmb_viirs_100m_2016.tif", pop_pix)

# Malaria incidence in 2018 (Malaria Atlas Project)
pop_pix$malaria <-
  get_cov("data/Zambia/Covariates/202206_Global_Pf_Incidence_Rate_ZMB_2018.tiff", pop_pix)

pop_pix$urban_prob <-
  get_cov("data/Zambia/Covariates/zmb_urb_prob.tif", pop_pix)


saveRDS(pop_pix, 
        paste0("data/Zambia/zmb_X_pop.rds"))

X_pop <- readRDS(paste0("data/Zambia/zmb_X_pop.rds"))

# transform
X_pop$l1pntl <- log1p(X_pop$ntl)
X_pop$l1paccess <- log1p(X_pop$access)

cen_pop <- rast("data/Zambia/Population/zmb_ppp_2010.tif")
cen_pop_1km <- aggregate(cen_pop, fact = 10, fun = "sum", na.rm = T)
covs <- c("urban", "l1paccess", "l1pntl", "malaria", "urban_prob")

rast_list <- lapply(covs, function(x) {
  rasterize(X_pop, cen_pop_1km,
            field = x)
})
cov_rast <- do.call("c", rast_list)
names(cov_rast) <- covs
#cov_rast$urban <- as.factor(cov_rast$urban)

my_theme <- 
  theme(legend.position="bottom",
        legend.text=element_text(size=12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border=element_blank(),
        axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        legend.key.size = unit(1, "cm")) 

#### Plot covariates

urban_gg <- ggplot() +
  geom_spatraster(data = cov_rast, aes(fill = urban)) +
  scale_fill_whitebox_c(palette = "viridi") +
  labs(fill = "", title = "Urban") +
  theme_minimal(base_size = 22)  + my_theme
urban_prob_gg <- ggplot() +
  geom_spatraster(data = cov_rast, aes(fill = urban_prob)) +
  scale_fill_whitebox_c(palette = "viridi") +
  labs(fill = "", title = "Urban probability") +
  theme_minimal(base_size = 22)  + my_theme
access_gg <- ggplot() +
  geom_spatraster(data = cov_rast, aes(fill = l1paccess)) +
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

cov_gg <- (urban_gg + access_gg) / (ntl_gg + malaria_gg)
ggsave("results/figures/Zambia_covariate_map.png",
       cov_gg,
       width = 10, height = 11.4)
