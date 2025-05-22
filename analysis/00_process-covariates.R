# devtools::install_github("richardli/surveyPrev")
library(surveyPrev)
library(rdhs)
library(raster)
library(sf)
library(patchwork)
library(terra)
library(spdep)
#---------------------------------------------------------------------------#
# Load data
#---------------------------------------------------------------------------#
dat <- readRDS("data/Zambia/clean_DHS_data.rds")

dat1 <- data.frame(dat[, c("cluster", "hshold", "wt", "urban", "hiv")])
dat1 <- dat1[, 1:5]
colnames(dat1) <- c("cluster", "householdID", "weight", "strata", "value")

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
poly_adm2$admin2.name.full <- paste(poly_adm2$NAME_1, poly_adm2$NAME_2, sep = "_")

# geo <- getDHSgeo(country = "Zambia", year = 2018)
geo <- read_sf("data/Zambia/DHS/ZMGE71FL/ZMGE71FL.shp")
cluster.info <- clusterInfo(geo=geo, poly.adm1=poly_adm1, poly.adm2=poly_adm2, 
                            by.adm1 = "NAME_1",by.adm2 = "NAME_2")

#---------------------------------------------------------------------------#
# Saved object later to avoid redoing computation...skip this chunk 
#---------------------------------------------------------------------------#
if(FALSE){
  raster <- raster(paste0("data/Zambia/Population/zmb_f_", 15, "_2018.tif"))
  for(i in seq(20, 45, by = 5)){
    tmp <- raster(paste0("data/Zambia/Population/zmb_f_", i, "_2018.tif"))
    raster <- raster + tmp
  }
  raster_2010 <- raster("data/Zambia/Population/zmb_ppp_2010.tif")
  agg.pop1 <- aggPopulation(tiff = raster,
                            poly.adm = poly.adm1,
                            by.adm = "NAME_1",
                            fact = 10)
  agg.pop2 <- aggPopulation(tiff = raster,
                            poly.adm = poly.adm2,
                            by.adm = "NAME_2",
                            by.adm.upper = "NAME_1",
                            fact = 10)
  
  # Not sure if this is the right census population table?? 
  # Just guessing from the file name
  table <- read.csv("data/Zambia/DHS/zmb_2018_poppa.csv")[1:10, ]
  urban.frac <- data.frame(admin1 = table[, 1], frac = table[, 6]/100)
  frac <- surveyPrev::getUR(tiff.census = raster_2010,
                            tiff.survey = raster,
                            prop.census = urban.frac,
                            fact = 10,
                            poly.adm1 =poly.adm1,
                            poly.adm2 =poly.adm2,
                            varname1 = "NAME_1",
                            varname2 = "NAME_2")
  frac$pop_survey <- wrap(frac$pop_survey)
  frac$UR_surface <- wrap(frac$UR_surface)
  save(agg.pop1, agg.pop2, frac, file = "data/Zambia/pop_info.RData")
}
#---------------------------------------------------------------------------#
# Saved object later to avoid redoing computation...skip this chunk 
#---------------------------------------------------------------------------#
if(FALSE){
  load("data/Zambia/pop_info.RData")
  
  ## set up covariates
  access <- terra::rast('data/Zambia/Covariates/2015_accessibility_to_cities_v1.0.tif')
  malaria <- terra::rast('data/Zambia/Covariates/202206_Global_Pf_Incidence_Rate_ZMB_2018.tiff')
  ntl <- terra::rast('data/Zambia/Covariates/zmb_viirs_100m_2016.tif')
  cov_raster_stack <- list(malaria = malaria,
                           access = access, l1paccess = log1p(access), 
                           ntl = ntl, l1pntl = log1p(ntl))
  zambia_cov <- surveyPrev::getCovariate(tiffs = cov_raster_stack,
                                         tiff.population = unwrap(frac$pop_survey),
                                         UR.surface = unwrap(frac$UR_surface),
                                         cluster.info = cluster.info,
                                         poly.adm = poly_adm2,
                                         by.adm ='NAME_2',
                                         by.adm.upper ='NAME_1',
                                         na.rm=TRUE,
                                         standardize = TRUE,
                                         # additional sum by 1, population already summed by fact = 10 before
                                         fact = 1)
  cov.unit <- 
    zambia_cov$cluster.cov[,c("cluster", "access", "l1paccess", "malaria", "ntl", "l1pntl")]
  cov.pixel <- 
    zambia_cov$natl.grid[,c("admin2.name.full","admin1.name","strata","Population",
                            "access", "l1paccess", "malaria", "ntl", "l1pntl")]
  save(zambia_cov, cov.unit, cov.pixel, file = "data/Zambia/cov_info.RData")
}

