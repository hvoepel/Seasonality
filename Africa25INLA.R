####################################################################
## Script: Africa25INLA.R
## Purpose: Spatiotemporal Bayesian modelling for 25 African countries
##     with the Bayesian package R: INLA.  Models use formula 'frml'
##     from linear regression model 'lm3' and removes flagged 'outliers'
##     from data table prior to model run.  Posterior marginals of fixed
##     covariates from initial run 'inla1' are used as prior parameters
##     in the subsequent run 'inla2', each run within the same function.
##     A lookup table is provided for covariate names and descriptions.
## NOTE TO USER: Data tables and model object cannot be provided to the
##     code user as per restrictions on data usage by the data source.
## -----------------------------------------------------------------
## Input: Model list cotaining 'lm3' ModelsFR (list) & ModelsTO (list),
##     corresponding data tables Africa25L1fr (tbl) & Africa25L1to (tbl),
##     and list of country-level adjacency matrices AdjMat (list). 
## Output: Updated ModelsFR and ModelsTO containing 'inla1' and 'inla2'
## -----------------------------------------------------------------
## Author: Hal E Voepel
## Date Created: Wed 27 Oct 2021
## Last Updated: Wed 15 Dec 2021
## Copyright (c) Hal E Voepel
## Email: H.E.Voepel@soton.ac.uk
## -----------------------------------------------------------------
## Sections: 
##    0: FUNCTION BLOCK
##    1: INLA modelling over 25 African countries
## 
## -----------------------------------------------------------------
## References: 
##    [1] https://becarioprecario.bitbucket.io/inla-gitbook/index.html
##    [2] https://www.paulamoraga.com/book-geospatial/index.html
## 
####################################################################

# Setting working directory and cleaning environment
# setwd('~/Seasonality')
rm(list = ls())

# Loading packages
library(tidyverse)
library(sf)
library(INLA)
library(INLAutils)
library(brinla)

# Loading data 
load('ModelsLM3.RData') # Two lists of 'lm1', 'lm2', and 'lm3' models for 'FR' and 'TO'
load('ModelDataAfrica25L1.RData') # AdjMat, lookup table, and 'fr' and 'to' data sets
load('gadm_africa25_sf.RData') # Level 0,1,2 admin polygons of 25 Africa countries
gadmL1 <- gadm_africa25_1 # Only need L1
rm(list = ls(pattern = 'gadm_'))



##==================================================================
## BEGIN FUNCTION BLOCK 
##------------------------------------------------------------------

ModelINLA <- function(lstModel, tabData, adjMat, iso3, prec = 0.001, quantiles = c(0.025, 0.5, 0.975)) {
  # lstModel  = named list of linear models for 25 African countries
  # tabData   = tibble dataframe containing data for 25 African countries
  # adjMat    = named list of adjacency matrices for 25 African countries
  # iso3      = country code scalar for modelling single country in function
  # prec      = precision on fixed parameters for initial model
  # quantiles  = quantile values calculated for fixed posterior values
  
  # Building country specific data table dropping NA values
  Afr <- tabData %>% filter(GID_0 == iso3) %>% drop_na()
  
  # Removing high leverage, high influence outlier data
  outliers <- lstModel[[iso3]]$lm3$outliers
  if (length(outliers) > 0) {
    Afr <- Afr[-outliers,]
  }
  
  # Centering country specific covariates
  vars <- Afr %>% select(days_HOLS:max_GDPC) %>% names()
  Afr <- Afr %>% 
    mutate(across(all_of(vars), scale, scale = FALSE)) %>% 
    mutate(across(all_of(vars), as.vector)) %>% as.data.frame()
   
  # Selecting country specific adjacency matrix
  adjMat <- adjMat[[iso3]]
  
  # Constructing formula string
  frml <- lstModel[[iso3]]$lm3$frml # Retrieve 'frml' 
  frmlST <- "+  
  f(ID1,
      model = 'bym2', 
      graph = adjMat) +
  f(month2, 
      model = 'ar1') +
  f(ID2, 
      model = 'bym2', 
      graph = adjMat,
      group = month, 
      control.group = list(model = 'ar1'))"
  frmlINLA <- paste(frml, frmlST) %>% 
    str_replace_all("[\r\n]" , "") %>% 
    str_squish()
  
  # ****************************************
  # Running initial INLA model: inla1
  lstModel[[iso3]]$inla1$frml <- frmlINLA
  lstModel[[iso3]]$inla1$model <- inla(
    formula = as.formula(frmlINLA),
    data = Afr,
    family = 'gaussian',
    control.fixed = list(prec = prec, quantiles = quantiles),
    control.predictor = list(compute = TRUE),
    control.compute = list(dic = TRUE), 
    control.lincomb = list(verbose = FALSE)
  )
  
  # Using 'inla1' posteriors to update 'inla2' priors
  # =================================================
  m <- lstModel[[iso3]]$inla1$model # reassigning model to variable
  
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
  # Running updated INLA model: inla2
  lstModel[[iso3]]$inla2$frml <- frmlINLA
  lstModel[[iso3]]$inla2$model <- inla(
    formula = as.formula(frmlINLA),
    data = Afr,
    family = 'gaussian',
    control.fixed = fixed.priors,
    control.family = list(hyper = list(prec = list(prior = table.string))),
    control.predictor = list(compute = TRUE),
    control.compute = list(dic = TRUE), 
    control.lincomb = list(verbose = FALSE)
  )
  
  # Returning named list of models
  return(lstModel)
  
}

## END FUNCTION BLOCK
##==================================================================



####################################################################
## 1: INLA modelling over 25 African countries
##------------------------------------------------------------------
## Summary: 
## Var (Type): 
####################################################################

# Increasing memory limit (Windows only)
# memory.limit()
# memory.limit(size = 159800) # 10 x default value

# Loading data (reload on each run)
load('ModelsLM3.RData')

# Getting 'GID_0' list of 25 African countries and display in console
(iso3 <- sort(unique(Africa25L1fr$GID_0)))

# Calculate precision for fixed terms prior: either numeric or specified prior
prec <- 0.001 

# Setting quantiles for marginals to include 90%CI
quantiles <- c(0.025, 0.05, 0.5, 0.95, 0.975)

# Setting processing & lapse timer
prc <- proc.time()

# Running INLA model for each country
for (i in iso3) {
  
  # Resetting lapse timer
  lps <- proc.time()
  
  # Print process to console
  print(paste('processing INLA for',i,'...'))
  
  # Running INLA models
  # Returns two INLA models: initial model (inla1), refined model (inla2)
  ModelsFR <- ModelINLA(ModelsFR, Africa25L1fr, AdjMat, i, prec, quantiles)
  ModelsTO <- ModelINLA(ModelsTO, Africa25L1to, AdjMat, i, prec, quantiles)
  
  # Printing lapse time to console
  print(paste('... lapse time for', i, 'was'))
  print(proc.time() - lps)
  print('********************************')
  
}

# Printing total processing time
print('Total processing time was')
print(proc.time() - prc)
 


# Checking size and saving results
lobstr::obj_size(ModelsFR)/1024^3 # as GB
lobstr::obj_size(ModelsTO)/1024^3 # as GB
# Removed outliers = 'out', see above for 'Run<no>' labels
save(ModelsFR, ModelsTO, file = 'Models_Final.RData')





# *******************************************************
# manual model run
ModelsTMP <- ModelINLA(ModelsFR, Africa25L1fr, AdjMat, 'BFA', prec, quantiles)
summary(ModelsTMP$BFA$inla1$model)$fixed

summary(ModelsFR$BFA$lm3$model)
summary(ModelsFR$BFA$inla1$model)$fixed

## End manual model run
##------------------------------------------------------------------

