# Code for "The Two Cultures of Prevalence Mapping: Small Area Estimation and Model-Based Geostatistics"

## Introduction

All code is written in R and uses packages available on CRAN. The testing version of INLA needs to be downloaded from <https://www.r-inla.org/download-install> and installed.

The paper reviews methods for estimating prevalences at the subnational level. This repository contains the code for a case study involving estimation of HIV prevalence among women aged 15–49 in Zambia, based on household survey data from the Zambia 2018 Demographic Health Survey.

## HIV prevalence in Zambia, 2018

We compare and contrast small area estimation and model-based geostatistics approaches using data from the 2018 Zambia DHS to estimate HIV prevalence among women aged 15-49. In our example, we characterize variation across 10 provinces (Admin-1 areas) and 72 districts (Admin-2 areas) of Zambia, based on the 2018 Demographic Health Survey (DHS), which used stratified (by urban/rural status crossed with Admin-1 area) two-stage cluster sampling.

### Code

* analysis/functions.R: Provides functions for fitting models.
- analysis/00_process-data.R: Loads and processes DHS data and covariates.
- analysis/01_get_estimates.R: Main script. Runs all models. 
- analysis/02_admin1-cross-validation.R: Leave-one-out cross validation at Admin-1 level.
- analysis/03_admin2-cross-validation.R: Leave-one-out cross validation at Admin-2 level.
- analysis/04_figures.R: Creates figures and tables for manuscript.


### Data

#### DHS

-   Sampling frame information for 2018 DHS survey: <https://dhsprogram.com/publications/publication-fr361-dhs-final-reports.cfm>

-   Zambia 2018 DHS survey (NDHS2018) + displaced GPS coordinates: <https://dhsprogram.com/>

#### WorldPop

* Population 2010: https://hub.worldpop.org/geodata/summary?id=3963

* Population (15-49 years, female) 2018: https://hub.worldpop.org/geodata/summary?id=16429

* Night time lights: https://hub.worldpop.org/geodata/summary?id=18778

* License: https://creativecommons.org/licenses/by/4.0/

#### The Malaria Atlas Project

-   Accessibility to cities 2015: <https://malariaatlas.org/research-project/accessibility-to-cities/>

* Malaria incidence 2018: https://data.malariaatlas.org

-   License: <http://creativecommons.org/licenses/by/4.0/>

#### GADM

-   Shape files of Zambia (versions 3.6 and 4.1; admin0, admin1 and admin2): <https://gadm.org/download_country.html>

-   License: Academic use only. No commercial use or redistribution.

#### Urban/Rural covariate surface
 
* Produced by Yunhan Wu: Available upon request.

