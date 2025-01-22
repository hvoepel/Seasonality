####################################################################
## Script: Africa25LM.R
## Purpose: Country-level multivariate regressions including 2-way
##     interaction terms to determine best significant covariate 
##     set for use in the spatiotemporal Bayesian model, INLA [1,2].
##     Country-level models regress log10 transformed baseline monthly
##     changes in mobility (relative to previous January) for 25 countries.
##     Input covariate and baseline change in mobility data were cleaned
##     and tables formatted for modelling in previous code (not included).
##     Two data tables for different mobility: source (FR) & destination (TO)
## NOTE TO USER: Data tables and model object cannot be provided to the
##     code user as per restrictions on data usage by the data source.
## -----------------------------------------------------------------
## Input: Clean data table of covariates and baseline change in mobility
## Output: 'ModelFR' (list) and 'ModelTO' (list) each containing linear
##         model outputs and performance metrics: 'lm1', 'lm2', and 'lm3'
## -----------------------------------------------------------------
## Author: Hal E Voepel
## Date Created: Mon 8 Nov 2021
## Last Updated: Wed 30 Mar 2022
## Copyright (c) Hal E Voepel
## Email: H.E.Voepel@soton.ac.uk
## -----------------------------------------------------------------
## Sections: 
##    0: FUNCTION BLOCK
##    1: Stepwise regression to select model terms for 'lm1' and 'lm2'
##    2: Manual reduction of stepwise models as 'lm3' model formula 'frml'
##    3: Regress final model using reduced formula 'frml' from 'lm3'
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

library(MASS) # Support functions and datasets
library(tidyverse) # easy data manipulation
library(caret) # easy machine learning workflow
library(leaps) # computing stepwise regression
library(ggfortify) # regression diagnostics

# Loading data
load('ModelDataAfrica25L1.RData')
load('gadm_africa25_sf.RData')

# Retrieving spatial L1 admin polygons
adm1 <- gadm_africa25_1
rm(list = ls(pattern = 'gadm'))

# Getting list of 25 African country ISO3 codes
(iso3 <- Africa25L1fr %>% pull(GID_0) %>% unique())



##==================================================================
## BEGIN FUNCTION BLOCK 
##------------------------------------------------------------------

# Elements not in set function: %!in%
'%!in%' <- function(x,y)!('%in%'(x,y))

# Variance Inflation Factor function: VIF
VIF <- function(model) {
  tryCatch(
    expr = {
      VIF <- performance::check_collinearity(model)
    },
    error = function(e){} # Don't print error to console
  )
  return(VIF)
}

# Variance Decomposition Proportions function: VCP
VDP <- function(model) {
  tryCatch(
    expr = {
      VDP <- mctest::eigprop(model)
    },
    error = function(e){} # Don't print error to console
  )
  return(VDP)
}

# Adjusted R2 for test data
AdjR2 <- function(model, test.data) {
  # Calculate predicted values
  pred <- predict(model, test.data)
  # Get lengths obs and covariates
  n <- nrow(test.data)
  p <- str_subset(names(coef(model)), pattern = ':')
  p <- length(setdiff(names(coef(model)),p)[-1])
  # Return estimated adjusted R squared
  return(1-(1-R2(pred, test.data$Y))*(n - 1)/(n - p))
}


# Predicted residual sums of squares: PRESS
PRESS <- function(model) {
  # calculate the predictive residuals
  pr <- residuals(model)/(1-lm.influence(model)$hat)
  # calculate the PRESS
  PRESS <- sum(pr^2)
  return(PRESS)
}

# Predicted R-squared: PredR2
PredR2 <- function(model) {
  # Use anova() to get the sum of squares for the linear model
  lm.anova <- anova(model)
  # Calculate the total sum of squares
  tss <- sum(lm.anova$'Sum Sq')
  # Calculate the predictive R^2
  pred.r.squared <- 1-PRESS(model)/(tss)
  return(pred.r.squared)
}

# Model fit statistics: ModelFitStats
ModelFitStats <- function(model) {
  r.sqr <- summary(model)$r.squared
  adj.r.sqr <- summary(model)$adj.r.squared
  ratio.adjr2.to.r2 <- adj.r.sqr / r.sqr
  pre.r.sqr <- PredR2(model)
  ratio.pre.r.sqr.to.r2 <- pre.r.sqr / r.sqr
  press <- PRESS(model)
  return.df <- data.frame(
    "R-squared" = r.sqr,
    "Adj R-sqr" = adj.r.sqr,
    "Pred R-sqr" = pre.r.sqr,
    "Adj.R2-to-R2" = ratio.adjr2.to.r2,
    "Pred.R2-to-R2" = ratio.pre.r.sqr.to.r2,
    "PRESS" = press
  )
  return(round(return.df,3))
}

# Cleaning formula to exclude NA coefficients
FrmlClean <- function(regModel) {
  # Retrieving coefficient names that result in NA values
  naCoef <- names(coef(regModel)[is.na(coef(regModel))])
  frml <- paste(setdiff(names(coef(regModel)), naCoef)[-1], collapse = ' + ')
  return(paste('Y ~',frml))
}

# Sequential stepwise regressions (max 1 term for each 20 obs)
# Bidirectional stepwise followed by backward stepwise
# 'lm1' is two-stage bidirectional stepwise regression: 
# 1st stage is finds main effects, and 2nd stage uses 1st
# stage results to find main effects and 2-way interactions
# 'lm2' is a backwards stepwise regression of 'lm1' to
# reduce model to simplest, most significant form.
StepModel <- function(lstModel, frml, df, iso3) {
  # lstModel = list to store models and other outputs
  # frml = constructed formula of main effect model terms 
  # df = table of level L1 data for single country (iso3)
  # iso3 = country code that matches 'df' data table
  
  # Setting k-value for stepwise AIC threshold
  value <- tibble(
    p = .1*2^-(0:9),
    k = round(qchisq(p,1,lower.tail=FALSE),3)
  )
  
  ## Bidirectional stepwise: Use to detect potential interactions: lm1 
  ##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ## Stepwise in both directions, twice: (1) main effects, (2) mains plus interactions
  
  # Regressing initial model that includes all covariates
  initReg <- lm(formula = frml, data = df)
  
  # Removing NA coefficients and regressing again
  frml <- formula(FrmlClean(initReg))
  initReg <- lm(formula = frml, data = df)
  
  # Initializing threshold and index values
  nTerms <- Inf # number of terms in model
  nInter <- Inf # number of interaction terms
  idx <- 0      # indexing for k-values
  
  
  # ********************* MAIN EFFECTS STEPWISE *****************************
  stepBoth <- stepAIC(initReg, direction = 'both', trace = 0, k = value$k[1])
  
  # Iterating over k-values to not exceed max allowable number of model terms,
  # number of interaction terms does not exceed 4 (flexible), and as an emergency
  # halt, the index value, idx, does not exceed the number of k-values.
  while (nTerms > floor(nrow(df) / 20) & nInter > 4 & idx < nrow(value)) {
    
    # Iterating index
    idx <- idx + 1
    
    # ********************* INTERACTION TERMS STEPWISE *****************************
    stepIntr <- stepAIC(stepBoth, scope = . ~ .^2 , direction = 'both', 
                        trace = 0, k = value$k[idx]) # iterating over k-values
    
    # Retrieving number of terms
    nTerms <- length(names(coef(stepIntr))[-1])
    nInter <- length(str_subset(names(coef(stepIntr)), pattern = ':'))
    
  }
  
  # Storing bidirectional stepwise results as 'lm1' model
  lstModel[[iso3]]$lm1$kVal <- value[idx,]
  lstModel[[iso3]]$lm1$frml <- as.character(stepIntr$call)[2]
  lstModel[[iso3]]$lm1$terms <- names(coef(stepIntr))[-1]
  lstModel[[iso3]]$lm1$model <- stepIntr
  lstModel[[iso3]]$lm1$performance <- ModelFitStats(stepIntr)
  lstModel[[iso3]]$lm1$VIF <- VIF(stepIntr)
  lstModel[[iso3]]$lm1$VDP <- VDP(stepIntr)
  
  
  ## Backwards stepwise: Use to initially reduce interaction models: lm2 
  ##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  # Initializing values
  nTerms <- Inf # number of terms in model
  nInter <- Inf # number of interation terms
  ratio <- 1    # ratio of nTerms_lm2 / nTerms_lm1
  idx <- 0      # indexing for k-values
  
  # Initial model using final bidirectional formula
  frml <- as.formula(lstModel[[iso3]]$lm1$frml)
  initReg <- lm(formula = frml, data = df)
  
  # Iterate while ratio of number of terms in 'lm2' is greater than half 
  # the number of terms in 'lm1', and the index does not exceed k-value table
  while (ratio > 0.5 & idx < 10) {
    
    # Iterating index
    idx <- idx + 1
    
    # Stepwise in both directions
    stepBack <- stepAIC(initReg, direction = 'backward', trace = 0, k = value$k[idx])
    
    # Retrieving number of terms
    nTerms <- length(names(coef(stepBack))[-1])
    nInter <- length(str_subset(names(coef(stepBack)), pattern = ':'))
    
    # Calculating ratio of model term counts
    ratio <- nTerms / length(lstModel[[iso3]]$lm1$terms)
    
  }
  
  # Storing bidirectional stepwise results as 'lm2' model
  lstModel[[iso3]]$lm2$kVal <- value[idx,]
  lstModel[[iso3]]$lm2$frml <- as.character(stepBack$call)[2]
  lstModel[[iso3]]$lm2$terms <- names(coef(stepBack))[-1]
  lstModel[[iso3]]$lm2$model <- stepBack
  lstModel[[iso3]]$lm2$performance <- ModelFitStats(stepBack)
  lstModel[[iso3]]$lm2$VIF <- VIF(stepBack)
  lstModel[[iso3]]$lm2$VDP <- VDP(stepBack)
  
  # Returning list of models
  return(lstModel)
  
}


## END FUNCTION BLOCK
##==================================================================




####################################################################
## 1: Stepwise regression to select model terms
##------------------------------------------------------------------
## Summary: Building 25 linear models for each FR and TO models sets: 
##      'lm1' = bidirectional stepwise, 'lm2' = backwards stepwise
##      Bidirectional in two stages: (1) main effects, (2) interactions
##      NOTE: Covariates are mean centered prior to regressing models.
## Var (Type): ModelsFR (list) ModelsTO (list) containing 'lm1' & 'lm2'
####################################################################


# Initializing model lists
ModelsFR <- list()
ModelsTO <- list()


## Build list of models and their metrics
##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

for (i in iso3) { # 
  
  # Printing current model run
  print(paste0('Processing stepwise model for ', i, '...'))
  
  ## Building stepwise models for 25 'fr' baseline mobility models
  ##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  # Building country specific data table dropping NA values
  Afr <- Africa25L1fr %>% filter(GID_0 == i) %>% select(-c(ID1:n,domestic)) %>% drop_na()
  
  # Centering country specific covariates
  vars <- names(Afr)[-1]
  Afr <- Afr %>% 
    mutate(across(vars, scale, scale = FALSE)) %>% 
    mutate(across(vars, as.vector))
  
  # Construct formula and get initial regression
  frml <- formula(paste('Y ~', paste(names(Afr)[-1], collapse = ' + ')))
  
  # Building 'ModelsFR' for current 'iso3' country
  ModelsFR <- StepModel(ModelsFR, frml, Afr, i)
  
  ## Building stepwise models for 25 'to' baseline mobility models
  ##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  # Building country specific data table dropping NA values
  Afr <- Africa25L1to %>% filter(GID_0 == i) %>% select(-c(ID1:n,domestic)) %>% drop_na()
  
  # Centering country specific covariates
  vars <- names(Afr)[-1]
  Afr <- Afr %>% 
    mutate(across(vars, scale, scale = FALSE)) %>% 
    mutate(across(vars, as.vector))
  
  # Construct formula and get initial regression
  frml <- formula(paste('Y ~', paste(names(Afr)[-1], collapse = ' + ')))
  
  # Building 'ModelsTO' for current 'iso3' country
  ModelsTO <- StepModel(ModelsTO, frml, Afr, i)
  
}



####################################################################
## 2: Manual reduction of stepwise models to final model: 'lm3'
##------------------------------------------------------------------
## Summary: Check summaries and performance metrics for each 'lm2' and
##    possibly 'lm1' model (if 'lm2' has issues) to reduce to final
##    linear model 'lm3'.  Corrected issues during manual model 'lm3'
##    reduction include multicollinearity, overfitting, unacceptable
##    quality metrics, and model performance.  Covariates for 'lm3'
##    models are then used in INLA spatiotemporal Bayesian modelling.
##    The final 'lm3' formula is stored for a final linear regression.
## Var (Type): Updated ModelsFR and ModelsTO to include 'lm3' formula
####################################################################



# Tabulating max number of model terms using min 20 obs / term
obs <- 20
TermCount <- Africa25L1fr %>% 
  group_by(GID_0) %>% 
  summarise(n_fr = n(), terms_fr = floor(n_fr/obs)) %>% 
  ungroup() %>% 
  inner_join(
    Africa25L1to %>% 
      group_by(GID_0) %>% 
      summarise(n_to = n(), terms_to = floor(n_to/obs)) %>% 
      ungroup(),
    by = 'GID_0'
  )

# Term count for each 'lm1' and 'lm2' model 
# NOTE: Attaching 'fr' and 'to' separately
nTerms <- tibble(
  lm1_fr = integer(length = 25),
  lm2_fr = integer(length = 25),
  lm1_to = integer(length = 25),
  lm2_to = integer(length = 25)
)

for (k in 1:25) {
  nTerms$lm1_fr[k] <- length(ModelsFR[[k]]$lm1$terms)
  nTerms$lm2_fr[k] <- length(ModelsFR[[k]]$lm2$terms)
  nTerms$lm1_to[k] <- length(ModelsTO[[k]]$lm1$terms)
  nTerms$lm2_to[k] <- length(ModelsTO[[k]]$lm2$terms)
}
TermCount <- bind_cols(TermCount, nTerms)
TermCount <- TermCount %>% relocate(c(lm1_fr,lm2_fr), .before = 'n_to')


## Tabulating performance metrics 
##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
mod <- ModelsFR # Assigning to short name, change 'fr' to 'to' as needed

# Constructing table of interaction terms with highest occurance across models
twoWayTerms <- character()
for (i in iso3) { # Collecting interaction terms from all models
  if (length(str_subset(mod[[i]]$lm1$terms, pattern = ':')) > 0) {
    twoWayTerms <- c(twoWayTerms, str_subset(mod[[i]]$lm1$terms, pattern = ':'))
  }
} # Tabulating
twoWay <- tibble(
  term = names(table(twoWayTerms)),
  count = unname(table(twoWayTerms))
) %>% arrange(desc(count))

# Tabulating all model performance metrics
tabPerform <- mod$AGO$lm1$performance
for (i in iso3[-1]) {
  tabPerform <- bind_rows(tabPerform, mod[[i]]$lm1$performance)
}
tabPerform <- bind_cols(tibble(GID_0 = iso3), tabPerform)

# Tabulating all (p, k) pairs
tabValues <- mod$AGO$lm1$kVal
for (i in iso3[-1]) {
  tabValues <- bind_rows(tabValues, mod[[i]]$lm1$kVal)
}
tabValues <- bind_cols(tibble(GID_0 = iso3), tabValues)


# Printing all tables to console
TermCount %>% print(n = nrow(.))
twoWay %>% filter(count > 1) %>% print(n = nrow(.))
tabPerform %>% arrange(Pred.R2.to.R2) %>% print(n = nrow(.))
tabValues %>% print(n = nrow(.))



##------------------------------------------------------------------
## 2.1:  Checking individual terms
####################################################################

# Outputs for 'lm1' and 'lm2'
mod <- ModelsFR
for (i in iso3[1:5 + 4*5]) {
  cat('****************************************************************\n')
  # cat('VDP summary for linear model...\n')
  # print(mod[[i]]$lm2$VDP)
  # cat('\n')
  # cat('================================================================\n')
  # cat('VIF for summary for linear model...\n')
  # print(mod[[i]]$lm2$VIF)
  # cat('\n')
  cat('================================================================\n')
  cat(paste('Summary of 2-way interation regression for',i, '\n'))
  cat('********** BIDIRECTIONAL STEPWISE ****************\n')
  print(summary(mod[[i]]$lm1$model))
  cat('\n')
  cat('\n')
  cat('\n')
  cat('********** BACKWARD STEPWISE **********************\n')
  print(summary(mod[[i]]$lm2$model))
  cat('\n')
  cat('================================================================\n')
  cat('Performance metrics for lm1 model...\n')
  print(mod[[i]]$lm1$modelPerform)
  cat('\n')
  cat('Performance metrics for lm2 model...\n')
  print(mod[[i]]$lm2$modelPerform)
  cat('\n')
  cat('Reduced model terms...\n')
  cat('********** BIDIRECTIONAL STEPWISE ****************\n')
  print(mod[[i]]$lm1$frml)
  cat('********** BAKCWARD STEPWISE **********************\n')
  print(mod[[i]]$lm2$frml)
  cat(paste('************************ END:',i, '******************************\n'))
  cat('\n')
  cat('\n')
  cat('\n')
  cat('\n')
  cat('\n')
}

# Plot diagnosis of linear models
for (i in iso3[1:5]) {
  for (k in 1:6) { # building plots
    plot(mod[[i]]$lm2$model, which = k, main = i)
  }
}


# Build table of terms ratios and corresponding k-values
ratio <- numeric(25)
kValue <- numeric(25)
for (k in 1:25) {
  ratio[k] <- round(length(mod[[iso3[k]]]$lm2$terms) / 
                      length(mod[[iso3[k]]]$lm1$terms), digits = 3)
  kValue[k] <- modelsINTRto[[iso3[k]]]$lm2$kVal$k
}
tabRatios <- tibble(
  ratio = ratio,
  kValue = kValue
) %>% mutate(GID_0 = iso3, .before = 'ratio') %>% 
  print(n = nrow(.))



## Manually select country-specific data, center, and drop outliers
##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## NOTE: Must change between 'fr' and 'to' in both 'Africa25..' and 'Models..'
## NOTE: Code used interactively in Console, not intended for sequential execution.

# Country specific data table
k <- 3
Afr <- Africa25L1fr %>% filter(GID_0 == iso3[k]) %>% select(-c(ID1:n,domestic)) %>% drop_na()

# Dropping outliers
(outliers <- ModelsFR[[iso3[k]]]$lm3$outliers)
if (length(outliers) > 0) {
  Afr <- Afr[-outliers,]
}

# Centering covariates
vars <- names(Afr)[-1]
Afr <- Afr %>% 
  mutate(across(vars, scale, scale = FALSE)) %>% 
  mutate(across(vars, as.vector))
# summary(Afr)

# Retriving lm3 and setup lm4 for reduction
lm2 <- lm(ModelsFR[[iso3[k]]]$lm2$frml, data = Afr)
lm3 <- update(lm2, . ~ .)
summary(lm3)

# Checking single covariate models
covTGO <- numeric()
idx <- 1
for (t in names(coef(lm3))[-1]) {
  frml <- formula(paste('Y ~',t))
  lm1 <- lm(frml, Afr)
  if (summary(lm1)$coef[8] < 0.1) {
    print(summary(lm1))
    covTGO[idx] <- t
    idx <- idx + 1
  }
}


## End manual selection
##------------------------------------------------------------------





####################################################################
## 3: Regress final model using reduced formula 'frml' from 'lm3'
##------------------------------------------------------------------
## Summary: Once final formula for 'lm3' is determined, use it to 
##     regress the final model including all quality and performance
##     metrics.  Check final models before moving on to INLA modelling.
## Var (Type): Updating ModelFR and ModelTO with 'lm3' run and metrics
####################################################################


# Model <- ModelsTO

for (i in iso3) {
  
  # Country specific data table
  Afr <- Africa25L1to %>% filter(GID_0 == i) %>% select(-c(ID1:n,domestic)) %>% drop_na()
  
  # Regress model for ith country removing unwanted outliers
  if (length(Model[[i]]$lm3$outliers) > 0) {
    out <- Model[[i]]$lm3$outliers
    frml <- as.formula(Model[[i]]$lm3$frml)
    lm3 <- lm(frml, Afr[-out,])
  } else {
    frml <- as.formula(Model[[i]]$lm3$frml)
    lm3 <- lm(frml, Afr)
  }
  
  # Centering country specific covariates
  vars <- names(Afr)[-1]
  Afr <- Afr %>% 
    mutate(across(vars, scale, scale = FALSE)) %>% 
    mutate(across(vars, as.vector))
  
  # Storing model, terms, and metrics
  Model[[i]]$lm3$terms <- names(coef(lm3))[-1]
  Model[[i]]$lm3$model <- lm3
  Model[[i]]$lm3$performance <- ModelFitStats(lm3)
  Model[[i]]$lm3$VIF <- VIF(lm3)
  
}

# # Reassign 'Model' and remove when finished
# ModelsTO <- Model
# rm(Model)


# Checking performance
# NOTE: change 'FR' to 'TO' and reassign afterwards
df <- data.frame(matrix(ncol = 6, nrow = 0))
names(df) <- names(ModelsTO$AGO$lm3$performance)
for (i in iso3) {
  df <- bind_rows(df, ModelsTO[[i]]$lm3$performance)
}
df <- bind_cols(data.frame(GID_0 = iso3), df)

# Reassigning
# perfFR <- df
# perfTO <- df

perfFR %>% arrange(Pred.R2.to.R2)
perfTO %>% arrange(Pred.R2.to.R2)

# Checking size and saving results
lobstr::obj_size(ModelsFR)/1024^3 # as GB
lobstr::obj_size(ModelsTO)/1024^3 # as GB
# save(ModelsFR, ModelsTO, file = 'Models.RData')


# Outputs for 'lm3'
mod <- modelsINTRfr
for (i in iso3) {
  cat('****************************************************************\n')
  # cat('VDP summary for linear model...\n')
  # print(mod[[i]]$lm3$VDP)
  # cat('\n')
  # cat('================================================================\n')
  # cat('VIF for summary for linear model...\n')
  # print(mod[[i]]$lm3$VIF)
  # cat('\n')
  cat('================================================================\n')
  cat(paste('Summary of 2-way interation regression for',i, '\n'))
  print(summary(mod[[i]]$lm3$model))
  cat('\n')
  cat('================================================================\n')
  cat('Performance metrics for linear model...\n')
  print(mod[[i]]$lm3$modelPerform)
  cat(paste('************************ END:',i, '******************************\n'))
  cat('\n')
  cat('\n')
  cat('\n')
  cat('\n')
  cat('\n')
}
