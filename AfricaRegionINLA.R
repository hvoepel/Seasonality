## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
## Script: Africa25INLA.R ------------------------------------------
## Purpose: Spatiotemporal Bayesian modelling using 'lm3' covariates
##     as fixed covariate input using the INLA modelling package[1]
##     and based on the African case study using INLA [2].  Coding here
##     is similar to 'Africa25INLA.R' script except modelled over regions.
##     A lookup table is provided for covariate names and descriptions.
## NOTE TO USER: Data tables and model object cannot be provided to the
##     code user as per restrictions on data usage by the data source.
##    
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
## Input: Data tables: AfricaRegionL1fr & AfricaRegionL1to
##        Model lists: RegionFR & RegionTO
## Output: Update RegionFR and RegionTO with 'inla1' and 'inla2' outputs
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
## Author: Hal E Voepel
## Date Created: Wed 15 Dec 2021
## Last Updated: Mon 28 Mar 2022
## Email: H.E.Voepel@soton.ac.uk
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
## Sections: 
##    0: Combining GADM tables as master tables and FUNCTION BLOCK
##    1: INLA modelling over 25 African countries by Region
## 
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
## References: 
##    [1] https://becarioprecario.bitbucket.io/inla-gitbook/index.html
##    [2] https://www.nature.com/articles/s41598-021-94683-7
## 
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# Setting working directory and cleaning environment
setwd('<PATH TO WORKSPACE FOLDER>')
rm(list = ls())


# Setting working directory and cleaning environment
# setwd('~/Seasonality') 
rm(list = ls())

# Loading packages
library(tidyverse)
library(sf)
library(spdep)
library(INLA)
# library(INLAutils)
# library(brinla)

# Loading data and model lists
load('ModelsRegional.RData')   
# "RegionFR"  "RegionTO"  "lookup"
load('ModelDataAfricaRegionL1.RData')
# "AfricaRegionL1fr" "AfricaRegionL1to"
# NOTE: Data are already ordered on GID_0, GID_1, month2
load('gadm_africa25_sf.RData')
# "gadm_africa25_0" "gadm_africa25_1" "gadm_africa25_2"
load('gadm_africaOther_sf.RData')
# "gadm_africaOther_0" "gadm_africaOther_1" "gadm_africaOther_2"

## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## 0: Combining GADM tables as master tables for parsing -----------
## *****************************************************************
## Summary: row bind two African GADM tables of disparate countries 
##    arrange level 0 (adm0master) on GID_0 and level 1 (adm1master) 
##    on GID_0, GID_1. Cast mixed geometry to MULTIPOLYGON geometry.
##
## Var (Type): adm0master (sf, tbl) adm1master (sf, tbl) 
##    
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Combining L0 tables and L1 tables and order 
adm0master <- bind_rows(gadm_africa25_0, gadm_africaOther_0) %>% arrange(GID_0)
adm1master <- bind_rows(gadm_africa25_1, gadm_africaOther_1) %>% arrange(GID_0, GID_1)
rm(list = ls(pattern = 'gadm')) # Cleaning up

# Casting to common geometry
st_geometry_type(adm0master)
adm0master <- st_cast(adm0master, 'MULTIPOLYGON') 
st_geometry_type(adm1master)
adm1master <- st_cast(adm1master, 'MULTIPOLYGON') 


## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
## *********************FUNCTION BLOCK******************************
## 
## SetAdminID: description 
## SetTableID: 
## SetAdjMat: 
## ModelINLA: 
## 
## =====================BEGIN FUNCTION BLOCK========================

# Sets IDs by level for specific data table: SetAdminID
SetAdminID <- function(tabData, level = 0) {
  
  if (level == 0) {
    L0 <- tabData %>% pull(GID_0) %>% unique()
    adm0master %>% filter(GID_0 %in% L0) %>% 
      arrange(GID_0) %>% 
      mutate(                 # ID.0.x = indexed on GID_0
        ID.0.0 = 1:nrow(.),
        ID.0.1 = ID.0.0,
        .before = 'GID_0'
      ) %>% return()
  } else {
    L1 <- tabData %>% pull(GID_1) %>% unique()
    adm1master %>% filter(GID_1 %in% L1) %>% 
      jamba::mixedSortDF('GID_1') %>%
      group_by(GID_0) %>%
      mutate(                 # ID.0.x = indexed on GID_1 grouped on GID_0
        ID.1.0 = seq(n()),
        ID.1.1 = ID.1.0,
        .before = 'GID_0'
      ) %>% ungroup() %>%
      mutate(                 # ID.0.x = indexed on GID_1 over entire Region
        ID.2.0 = seq(n()),
        ID.2.1 = ID.2.0,
        .before = 'GID_0'
      ) %>% return()
  }
  
}


# Updates data table with L0 and L1 admin IDs: SetTableID
SetTableID <- function(tabData, admLvl0, admLvl1) {
  
  # Setting Lvl1 IDs
  tabData <- admLvl1 %>% 
    st_drop_geometry() %>% 
    select(ID.1.0:GID_0, GID_1) %>%
    right_join(tabData, by = c('GID_0','GID_1'))
  
  # Setting Lvl0 IDs
  tabData <- admLvl0 %>% 
    st_drop_geometry() %>% 
    select(ID.0.0:GID_0) %>%
    right_join(tabData, by = 'GID_0') %>% 
    relocate(GID_0, .before = 'NAME_0') %>% 
    relocate(GID_1, .before = 'NAME_1') 
  
}


# Constructs adjacency matrix for admin sf: SetAdjMat
SetAdjMat <- function(admin) {
  
  # Convert to sp object, and process to adj.mat
  admin.sp <- as(st_geometry(admin),"Spatial") 
  admin.adj <- poly2nb(admin.sp)           
  adj.mat <- as(nb2mat(admin.adj, zero.policy = TRUE), "Matrix")
  return(adj.mat)
  
}


# INLA model fitting in two passes: ModelINLA
#   First pass initial fitting using uninformed priors: inla1
#   Second pass uses covariate posterior parameters as priors: inla2
ModelINLA <- function(lstModel, tabData, region, prec = 0.001, quantiles = c(0.025, 0.5, 0.975), pred = FALSE) {
  # lstModel  = named list of linear models for 25 African countries
  # tabData   = tibble dataframe containing data for 25 African countries
  # region    = region for modelling groups of countries in function
  # prec      = precision on fixed parameters for initial model
  # quantiles = quantile values calculated for fixed posterior values
  # pred      = whether to include prediction covariates (TRUE) or not (FALSE)
  
  # Building region specific data table
  # WARNING: Do not drop NA values; all predict covariates will be removed
  if (pred) {
    Afr <- tabData %>% filter(Region == region)         # For 'inla3' and 'inla4'
  }else {
    Afr <- tabData %>% filter(Region == region, n > 0)  # For 'inla1' and 'inla2'
  }
  
  # Setting 'domestic' as binary factor (this starts for Run5)
  # NOTE: 'domestic' cannot be used in predicting 
  Afr <- Afr %>% mutate(domestic = factor(domestic))
  
  # Removing high leverage, high influence outliers determine in previous run
  outliers <- lstModel[[region]]$lm3$outliers
  if (length(outliers) > 0) {
    Afr <- Afr[-outliers,]
  }
  
  # Removing problematic islands
  if (region == 'Eastern') {
    Afr <- Afr %>% filter(GID_0 != 'MDG')
  }
  
  # Setting IDs for admin levels to correspond with selected data table
  adm0 <- SetAdminID(Afr, level = 0)
  adm1 <- SetAdminID(Afr, level = 1)
  
  # Setting IDs for data tables from previously ID'd admin sets
  Afr <- SetTableID(Afr, adm0, adm1)
  
  # Centering region specific covariates
  vars <- Afr %>% select(days_HOLS:max_GDPC) %>% names()
  Afr <- Afr %>% 
    mutate(across(all_of(vars), scale, scale = TRUE)) %>% 
    mutate(across(all_of(vars), as.vector)) %>% as.data.frame()
  
  # Selecting region and country specific adjacency matrix
  adjMat0 <- SetAdjMat(adm0)
  adjMat1 <- SetAdjMat(adm1)
  
  # Constructing formula string
  frml <- lstModel[[region]]$lm3$frml # Retrieve 'frml'
  # f(GID_0) + 
  # f(domestic,
  #     model = 'iid') +
  # f(year, 
  #   model = 'ar1') + 
  # f(ID.0.1,
  #   model = 'bym2',
  #   graph = adjMat0,
  #   group = GID_1) + 
  # LEVELS ON ID (x = copy index)
  # ID.0.x = indexed on GID_0 (use with adjMat0)
  # ID.1.x = indexed on GID_1 grouped on GID_0 (use with adjMat1)
  # ID.2.x = indexed on GID_1 over entire Region (use with adjMat1)
  frmlST <- "+ 
  f(ID.0.0,
      model = 'bym2',
      graph = adjMat0) + 
  f(ID.1.0,
      model = 'bym2', 
      graph = adjMat1) + 
  f(ID.2.0,
      model = 'bym2', 
      graph = adjMat1) + 
  f(month2, 
      model = 'ar1') + 
  f(ID.2.1, 
      model = 'bym2', 
      graph = adjMat1,
      group = month2, 
      control.group = list(model = 'ar1'))"
  frmlINLA <- paste(frml, frmlST) %>%
    str_replace_all("[\r\n]" , "") %>%
    str_squish()
  
  # '+ domestic', copy this between 'frml' and 'frmlST' above as needed
  
  # ****************************************
  # Running initial INLA model: inla1 (inla3 for final model)
  cat('...fitting inla1 model...\n')
  lstModel[[region]]$inla1$frml <- frmlINLA
  lstModel[[region]]$inla1$model <- inla(
    formula = as.formula(frmlINLA),
    data = Afr,
    family = 'gaussian',
    control.fixed = list(prec = prec, quantiles = quantiles),
    control.predictor = list(compute = TRUE),
    control.compute = list(dic = TRUE, cpo = TRUE, waic = TRUE, 
                           return.marginals.predictor = TRUE), 
    control.inla = list(strategy = 'gaussian'),
    control.lincomb = list(verbose = FALSE)
  )
  
  # Using 'inla1' posteriors to update 'inla2' priors
  # =================================================
  m <- lstModel[[region]]$inla1$model # reassigning model to variable
  
  # Utility functions to extract posterior moments
  posterior.mean <- function(m, name) m$summary.fixed[name,'mean']
  posterior.prec <- function(m, name) 1/m$summary.fixed[name,'sd']^2
  
  # Preparing priors for 'inla2' fixed effects: fixed.priors
  fixed.priors <- list()
  fixed.priors$mean = list()
  fixed.priors$prec = list()
  fixed.priors$quantiles = quantiles
  for (name in m$names.fixed) {
    fixed.priors$mean[[name]] <- posterior.mean(m, name)
    fixed.priors$prec[[name]] <- posterior.prec(m, name)
  }
  
  # Preparing priors for 'inla2' error term: table.string
  imh <- m$internal.marginals.hyperpar[[1]] # NOTE: Must use internal marginals!
  table.string <- paste(c('table:', imh[,'x'], imh[,'y']), collapse = ' ')
  
  
  # ****************************************
  # Running updated INLA model: inla2 (inla4 for final model)
  cat('...fitting inla2 model...\n')
  lstModel[[region]]$inla2$frml <- frmlINLA
  lstModel[[region]]$inla2$model <- inla(
    formula = as.formula(frmlINLA),
    data = Afr,
    family = 'gaussian',
    control.fixed = fixed.priors,
    control.family = list(hyper = list(prec = list(prior = table.string))),
    control.predictor = list(compute = TRUE),
    control.compute = list(dic = TRUE, cpo = TRUE, waic = TRUE, 
                           return.marginals.predictor = TRUE), 
    control.inla = list(strategy = 'gaussian'),
    control.lincomb = list(verbose = FALSE)
  )
  
  # Returning named list of models
  return(lstModel)
  
}

## =====================END FUNCTION BLOCK==========================  


## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## 1: INLA modelling over African countries by Region --------------
## *****************************************************************
## Summary: INLA modelling in two passes: 'inla1' is initial run
##    with uniformed priors, 'inla2' uses posteriors from 'inla1' as
##    prior information for building 'inla2' models
## Var (Type): 
##    
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# # Increasing memory limit (Windows only)
# memory.limit()
# memory.limit(size = 159800) # 10 x default value

# Loading data (reload on each run)
load('ModelsRegional.RData')
# load('ModelsRegionalFinal.RData')

# Getting 'Region' list of 25 African countries and display in console
(region <- sort(unique(AfricaRegionL1fr$Region)))

# Calculate precision for fixed terms prior: either numeric or specified prior
prec <- 0.001 # Run1 = 0.001, Run2 = 0.01, Run3 = 0.1, Run4 = 0.5, Run5 = 0.8
# NOTE: Run4 and Run5 likely creating biased estimated of fixed posteriors

# Setting quantiles for marginals
quantiles <- c(0.025, 0.05, 0.5, 0.95, 0.975)

# Setting processing & lapse timer
prc <- proc.time()

# Running INLA model for each country
for (r in region[4]) {
  
  # Resetting lapse timer
  lps <- proc.time()
  
  # Print process to console
  cat(paste('processing INLA for',r,'based on lm3 results...\n'))
  
  # INLA(lstModel, tabData, region, prec = 0.001, quantiles = c(0.025, 0.5, 0.975), pred = FALSE)
  # lstModel  = named list of linear models for 25 African countries
  # tabData   = tibble dataframe containing data for 25 African countries
  # region    = region for modelling groups of countries in function
  # prec      = precision on fixed parameters for initial model
  # quantiles = quantile values calculated for fixed posterior values
  # pred      = whether to include prediction covariates (TRUE) or not (FALSE)
  
  # Running INLA models
  # Returns two INLA models: initial model (inla1), refined model (inla2)
  
  cat('running INLA source (FR) model...\n')
  RegionFR <- ModelINLA(RegionFR, AfricaRegionL1fr, r, prec, quantiles, pred = TRUE)
  cat('running INLA destination (TO) model...\n')
  RegionTO <- ModelINLA(RegionTO, AfricaRegionL1to, r, prec, quantiles, pred = TRUE)
  
  # Printing lapse time to console
  cat(paste('... lapse time for', r, 'was\n'))
  print(proc.time() - lps)
  cat('********************************\n')
  
}

# Printing total processing time
cat('Total processing time was\n')
print(proc.time() - prc)





##------------------------------------------------------------------
## 1.1: Checking data and saving to image
####################################################################

# Checking size and saving results
lobstr::obj_size(RegionFR)/1024^2 # as MB
lobstr::obj_size(RegionTO)/1024^2 # as MB

# Removed outliers = 'out', see above for 'Run<no>' labels
# save(RegionFR, RegionTO, file = 'Models_Regional_Predict2.RData')







