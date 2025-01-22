####################################################################
## Script: AfricaRegionLM.R
## Purpose: Multivariate linear regression to determine covariate set
##     for regional-level models.  Modelling steps noted in sections.
##     Coding is similar to 'Africa25LM.R' except modelled over regions.
## NOTE TO USER: Data tables and model object cannot be provided to the
##     code user as per restrictions on data usage by the data source.
## -----------------------------------------------------------------
## Input: Covariate data table constructed in 'DataMungingRegionL1.R'
## Output: RegionFR (list) RegionTO (list)
## -----------------------------------------------------------------
## Author: Hal E Voepel
## Date Created: Mon 8 Nov 2021
## Last Updated: Wed 30 Mar 2022
## Copyright (c) Hal E Voepel
## Email: H.E.Voepel@soton.ac.uk
## -----------------------------------------------------------------
## Sections: 
##    0: FUNCTION BLOCK
##    1: Stepwise selection of model terms
##    2: Evaluating model performance and analysis for issues
##    3: Manual reduction of 'lm3' models checking for issues
## 
## -----------------------------------------------------------------
## References: 
##    [1] https://becarioprecario.bitbucket.io/inla-gitbook/index.html
##    [2] https://www.nature.com/articles/s41598-021-94683-7
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
load('ModelDataAfricaRegionL1.RData')
# "AfricaRegionL1fr" "AfricaRegionL1to"

# Getting list of 25 African country ISO3 codes
(iso3 <- AfricaRegionL1fr %>% pull(GID_0) %>% unique() %>% sort())


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
  VDP <- 'NA'
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

# Checking test data against model
TestModel <- function(model, test.data) {
  pred.data <- predict(model, test.data)
  tab <- tibble(
    R2 = caret::R2(pred.data, test.data$Y),
    RMSE = caret::RMSE(pred.data, test.data$Y),
    MAE = caret::MAE(pred.data, test.data$Y),
    PER = RMSE/mean(test.data$Y)
  )
  return(round(tab,3))
}

# Multiple stepwise regressions (assume 1:20 rule)
# Bidirectional stepwise followed by backward stepwise
StepModel <- function(lstModel, frml, df, idx, region) {
  # lstModel = list to store model and other outputs
  # frml = formula of main effect terms from 'df'
  # df = tibble containing 'region' level data
  # region = Africa region that matches 'df' data
  # idx = indices of testing data set
  
  cat(paste('Processing models for the', region, 'region of Africa\n'))
  
  # Partition dataframe 'df' into training/testing sets via 'idx'
  # dfTrain <- df[idx,]
  # dfTest <- df[-idx,]
  
  ## Bidirectional stepwise: Use to detect potential interactions 
  ##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ## Stepwise in both directions, twice: (1) main effects, (2) interactions
  
  # Regressing initial model that includes all covariates
  initReg <- lm(formula = frml, data = df)
  
  # Removing NA coefficients and regressing again
  frml <- formula(FrmlClean(initReg))
  initReg <- lm(formula = frml, data = df)
  
  # ********************* MAIN EFFECTS STEPWISE *****************************
  cat('Performing bidirectional stepwise for main effects...\n')
  stepBoth <- stepAIC(initReg, direction = 'both', trace = 0)
  print(summary(stepBoth))
  
  cat('\n\n Performing backward stepwise for main effects...\n')
  stepBack <- stepAIC(stepBoth, direction = 'backward', trace = 0, k = 5)
  print(summary(stepBack))
  
  # ********************* INTERACTION TERMS STEPWISE *****************************
  cat('\n\n Performing bidirectional stepwise for interaction effects...\n')
  stepIntr <- stepAIC(stepBack, scope = . ~ .^2 , direction = 'both', trace = 0, k = 5)
  print(summary(stepIntr))
  
  # Storing bidirectional stepwise results as 'lm1' model
  cat('\n\n Storing model and quality check metrics as lm1\n')
  lstModel[[region]]$lm1$frml <- as.character(stepIntr$call)[2]
  lstModel[[region]]$lm1$terms <- names(coef(stepIntr))[-1]
  lstModel[[region]]$lm1$model <- stepIntr
  lstModel[[region]]$lm1$test <- TestModel(stepIntr, df)
  lstModel[[region]]$lm1$perf <- ModelFitStats(stepIntr)
  lstModel[[region]]$lm1$VIF <- VIF(stepIntr)
  lstModel[[region]]$lm1$VDP <- VDP(stepIntr)
  # lstModel[[region]]$lm1$Idx <- as.vector(idx)
  
  
  ## Backwards stepwise: Use to reduce interaction models 
  ##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  # Initial model using final 'lm1' formula
  frml <- as.formula(lstModel[[region]]$lm1$frml)
  initReg <- lm(formula = frml, data = df)
  
  # Backward stepwise regression of 'lm1' model
  cat('Performing backwards stepwise of lm1 model...\n')
  stepBack <- stepAIC(initReg, direction = 'backward', trace = 0, k = 15)
  print(summary(stepBack))
  
  # Building table of metrics to check model against test data
  # pred <- predict(stepBack, df)
  
  # Storing backward stepwise results as 'lm2' model
  cat('\n\n Storing model and quality check metrics as lm2\n')
  lstModel[[region]]$lm2$frml <- as.character(stepBack$call)[2]
  lstModel[[region]]$lm2$terms <- names(coef(stepBack))[-1]
  lstModel[[region]]$lm2$model <- stepBack
  lstModel[[region]]$lm2$test <- TestModel(stepBack, df)
  # lstModel[[region]]$lm2$pred <- pred
  lstModel[[region]]$lm2$perf <- ModelFitStats(stepBack)
  lstModel[[region]]$lm2$VIF <- VIF(stepBack)
  lstModel[[region]]$lm2$VDP <- VDP(stepBack)
  # lstModel[[region]]$lm2$Idx <- as.vector(idx)
  
  cat('\n\n\n')
  
  # Returning list of models
  return(lstModel)
  
}


## END FUNCTION BLOCK
##==================================================================



####################################################################
## 1: Stepwise selection of model terms
##------------------------------------------------------------------
## Summary: Two bidirectional stepwise regressions, first for main
##     effects, the second for interaction terms.  Then a backward
##     stepwise regression is performed to simplify the covariate set.
## Var (Type): RegionFR (list)  RegionTO (list) with 'lm1' and 'lm2'
####################################################################


# Checking counts by region
Africa25L1fr %>% 
  group_by(Region) %>% 
  summarise(n = n())


## Build list of models and their metrics
##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Building single linear model for each models set: RegionTO, RegionFR
## 'lm1' = bidirectional stepwise, 'lm2' = backwards stepwise
## Bidirectional is in two stages: (1) main effects, (2) interactions
## NOTE: Covariates are mean centered prior to regressing models


## Building regional stepwise mobility models for 'fr' and 'to'
##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Getting list of Africa regions
(regions <- levels(Africa25L1fr$Region))

# Initializing model lists 
RegionFR <- list()
RegionTO <- list()

for (r in regions) {
  
  # Building country specific data tables dropping NA values
  AfrFR <- Africa25L1fr %>% filter(Region == r) %>% select(-c(ID1:n,domestic)) %>% drop_na()
  AfrTO <- Africa25L1to %>% filter(Region == r) %>% select(-c(ID1:n,domestic)) %>% drop_na()
  
  # Build models on both centering and non-centering (comment out for latter)
  varsFR <- names(AfrFR)[-1]
  varsTO <- names(AfrTO)[-1]
  AfrFR <- AfrFR %>% # Centering covariates only
    mutate(across(all_of(varsFR), scale, scale = TRUE)) %>%
    mutate(across(all_of(varsFR), as.vector))
  AfrTO <- AfrTO %>% # Centering covariates only
    mutate(across(all_of(varsTO), scale, scale = TRUE)) %>%
    mutate(across(all_of(varsTO), as.vector))
  
  # Partitioning (with/without centering) into train/test data sets
  idxFR <- createDataPartition(AfrFR$Y, p = 0.7, list = FALSE)
  idxTO <- createDataPartition(AfrTO$Y, p = 0.7, list = FALSE)
  
  # Construct formula and get initial regression
  (frmlFR <- formula(paste('Y ~', paste(names(AfrFR)[-1], collapse = ' + '))))
  (frmlTO <- formula(paste('Y ~', paste(names(AfrTO)[-1], collapse = ' + '))))
  
  # Building predictive models for 'fr' and 'to' Africa mobility
  RegionFR <- StepModel(RegionFR, frmlFR, AfrFR, idxFR, r)
  RegionTO <- StepModel(RegionTO, frmlTO, AfrTO, idxTO, r)
  
}



####################################################################
## 2: Evaluating model performance and analysis for issues
##------------------------------------------------------------------
## Summary: Construct model performance tables to use for manually
##   reducing to final 'lm3' model, which is used in INLA modelling
## Var (Type): N/A
####################################################################



# Retrieving model performance metrics
modPerformanceFR <- tibble(
  Region = character(),
  Model = character(),
  R.squared = numeric(),
  Adj.R.sqr = numeric(),
  Pred.R.sqr = numeric(),
  Adj.R2.to.R2 = numeric(),
  Pred.R2.to.R2 = numeric(),
  PRESS = numeric()
) 
modPerformanceTO <- modPerformanceFR

# Building separate model performance tables for 'FR' and 'TO'
for (r in regions) {
  for (m in c('lm1','lm2')) {
    tmp <- bind_cols(tibble(Region = r, Model = m),RegionFR[[r]][[m]]$perf,RegionFR[[r]][[m]]$test)
    modPerformanceFR <- bind_rows(modPerformanceFR, tmp)
    tmp <- bind_cols(tibble(Region = r, Model = m),RegionTO[[r]][[m]]$perf,RegionTO[[r]][[m]]$test)
    modPerformanceTO <- bind_rows(modPerformanceTO, tmp)
  }
}

# Combine all tables into single table
modPerformance <- bind_rows(
  bind_cols(
    tibble(
      Dir = rep('FR',10)
    ),
    modPerformanceFR
  ),
  bind_cols(
    tibble(
      Dir = rep('TO',10)
    ),
    modPerformanceTO
  )
) %>% arrange(Region, Dir, Model) #%>% 
  # mutate(Test.R2.to.R2 = round(R2 / R.squared, 3))

# Cleaning up
rm(modPerformanceFR, modPerformanceTO)
modPerformance




####################################################################
## 3: Manual reduction of 'lm3' models checking for issues
##------------------------------------------------------------------
## Summary: Similar to manual reduction of 'lm2' as 'lm3' that was
##    processed in 'Africa25LM.R'.  Coding is in snippets and not
##    intended to be run sequentially.
## Var (Type): RegionFR (list) and RegionTO (list) updated with 'lm3'
####################################################################



# Checking variability of GID_0 by region
# High variability in mean and variance
library(patchwork)

p <- list()
theme_set(
  theme_bw() +
    theme(legend.position = 'top')
)
for (r in regions) {
  
  p[[r]] <- Africa25L1fr %>% 
    filter(Region == r) %>% 
    ggplot(aes(x = Y, y = GID_0)) +
    geom_boxplot() +
    labs(title = r) +
    ylab('') + xlab('log baseline mobility')

}

# Arrange panel of plots
# (p$Northern / p$Central / p$Southern) | (p$Western / p$Eastern)

p$Western | p$Eastern | (p$Northern / p$Central / p$Southern)



# "Northern" "Central"  "Western"  "Eastern"  "Southern"
(r <- regions[1])

# Checking performance
modPerformance %>% 
  select(-c(Adj.R.sqr,Adj.R2.to.R2,PRESS,RMSE:PER)) %>%
  filter(Model == 'lm2', Region == r)

# Building country specific data tables dropping NA values
AfrFR <- Africa25L1fr %>% filter(Region == r) %>% select(-c(ID1:n,domestic)) %>% drop_na()
AfrTO <- Africa25L1to %>% filter(Region == r) %>% select(-c(ID1:n,domestic)) %>% drop_na()

# Build models on both centering and non-centering (comment out for latter)
varsFR <- names(AfrFR)[-1]
varsTO <- names(AfrTO)[-1]
AfrFR <- AfrFR %>% # Centering covariates only
  mutate(across(all_of(varsFR), scale, scale = TRUE)) %>% 
  mutate(across(all_of(varsFR), as.vector))
AfrTO <- AfrTO %>% # Centering covariates only
  mutate(across(all_of(varsTO), scale, scale = TRUE)) %>% 
  mutate(across(all_of(varsTO), as.vector))

# Checking both tables 
summary(AfrFR)
summary(AfrTO)




##------------------------------------------------------------------
## 3.1: Manual reduction of complex stepwise models 
####################################################################

# Tabulating max number of model terms using min 15 obs / term
obs <- 20
(TermCount <- Africa25L1fr %>% 
  group_by(Region) %>% 
  summarise(n_fr = n(), terms_fr = floor(n_fr/obs)) %>% 
  ungroup() %>% 
  inner_join(
    Africa25L1to %>% 
      group_by(Region) %>% 
      summarise(n_to = n(), terms_to = floor(n_to/obs)) %>% 
      ungroup(),
    by = 'Region'
  ))

# Term count for each 'lm1' and 'lm2' model 
# NOTE: Attaching 'fr' and 'to' separately
nTerms <- tibble(
  lm1_fr = integer(length = 5),
  lm2_fr = integer(length = 5),
  lm1_to = integer(length = 5),
  lm2_to = integer(length = 5)
)

for (k in 1:5) {
  nTerms$lm1_fr[k] <- length(RegionFR[[k]]$lm1$terms)
  nTerms$lm2_fr[k] <- length(RegionFR[[k]]$lm2$terms)
  nTerms$lm1_to[k] <- length(RegionTO[[k]]$lm1$terms)
  nTerms$lm2_to[k] <- length(RegionTO[[k]]$lm2$terms)
}
TermCount <- bind_cols(TermCount, nTerms)
TermCount <- TermCount %>% relocate(c(lm1_fr,lm2_fr), .before = 'n_to')
TermCount




## 3.1.1: Tabulating performance metrics 
##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
mod <- RegionFR # Assigning to short name, change 'fr' to 'to' as needed

# Constructing table of interaction terms with highest occurance across models
twoWayTerms <- character()
for (i in iso3) { # Collecting interaction terms from all models
  if (length(str_subset(mod[[i]]$lm1$terms, pattern = ':')) > 0) {
    twoWayTerms <- c(twoWayTerms, str_subset(mod[[i]]$lm1$terms, pattern = ':'))
  }
} 

# Tabulating all model performance metrics
tabPerform <- mod$Northern$lm1$perf
for (r in regions[-1]) {
  tabPerform <- bind_rows(tabPerform, mod[[i]]$lm1$perf)
}
tabPerform <- bind_cols(tibble(Region = regions), tabPerform)
tabPerform



# Printing all tables to console
TermCount %>% print(n = nrow(.))
tabPerform %>% arrange(Pred.R2.to.R2) %>% print(n = nrow(.))




## 3.1.2: Manual reduction of 'lm2' (or possibly 'lm1') models to 'lm3'
##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## "Northern" "Central"  "Western"  "Eastern"  "Southern"

library(tidyverse)

(r <- regions[3])

# Building country specific data tables dropping NA values
AfrFR <- Africa25L1fr %>% filter(Region == r) %>% select(-c(ID1:n,domestic)) %>% drop_na()
AfrTO <- Africa25L1to %>% filter(Region == r) %>% select(-c(ID1:n,domestic)) %>% drop_na()

# Build models on both centering and non-centering (comment out for latter)
varsFR <- names(AfrFR)[-1]
varsTO <- names(AfrTO)[-1]
AfrFR <- AfrFR %>% # Centering covariates only
  mutate(across(all_of(varsFR), scale, scale = TRUE)) %>%
  mutate(across(all_of(varsFR), as.vector))
AfrTO <- AfrTO %>% # Centering covariates only
  mutate(across(all_of(varsTO), scale, scale = TRUE)) %>%
  mutate(across(all_of(varsTO), as.vector))

# Manual reduction of FR...

summary(RegionFR[[r]]$lm3$model)
lm2 <- lm(RegionFR[[r]]$lm2$frml, AfrFR)
lm3 <- update(lm2, .~.)
RegionFR[[r]]$lm3$frml <- "Y ~ days_HSCH + mean_INFM * max_PRES + max_PRES * max_ACCS + mean_NTLS * max_PRES"
RegionFR[[r]]$lm3$frml <- RegionFR[[r]]$lm2$frml

# Manual reduction of TO...

summary(RegionTO[[r]]$lm2$model)
lm2 <- lm(RegionTO[[r]]$lm2$frml, AfrTO)
lm3 <- update(lm2, .~.)
RegionTO[[r]]$lm3$frml <- "Y ~ mean_EVAP + mean_NPRF + mean_SECF + mean_INFM * max_PRES"
RegionTO[[r]]$lm3$frml <- RegionTO[[r]]$lm2$frml

# Checking 'lm3' formulas
for (r in regions) {
  cat(paste('Formulas for lm3 at',r,'\n'))
  cat('RegionFR...\n')
  print(RegionFR[[r]]$lm3$frml)
  cat('RegionTO...\n')
  print(RegionTO[[r]]$lm3$frml)
  cat('\n\n')
}


# Fitting manually reduced 'lm3' models
for (r in regions) {
  
  cat(paste('Regressing lm3 models for the',r,'region...\n'))
  
  # Building country specific data tables dropping NA values
  AfrFR <- Africa25L1fr %>% filter(Region == r) %>% select(-c(ID1:n,domestic)) %>% drop_na()
  AfrTO <- Africa25L1to %>% filter(Region == r) %>% select(-c(ID1:n,domestic)) %>% drop_na()
  
  # Build models on both centering and non-centering (comment out for latter)
  varsFR <- names(AfrFR)[-1]
  varsTO <- names(AfrTO)[-1]
  AfrFR <- AfrFR %>% # Centering covariates only
    mutate(across(all_of(varsFR), scale, scale = TRUE)) %>%
    mutate(across(all_of(varsFR), as.vector))
  AfrTO <- AfrTO %>% # Centering covariates only
    mutate(across(all_of(varsTO), scale, scale = TRUE)) %>%
    mutate(across(all_of(varsTO), as.vector))
  
  # Making formulas
  frmlFR <- as.formula(RegionFR[[r]]$lm3$frml)
  frmlTO <- as.formula(RegionTO[[r]]$lm3$frml)
  
  # Regressing 'lm3' models
  lm3FR <- lm(frmlFR, AfrFR)
  lm3TO <- lm(frmlTO, AfrTO)
  
  # Storing 'FR' model results and metrics as 'lm3'
  cat('Storing FR model and quality check metrics as lm3\n')
  RegionFR[[r]]$lm3$terms <- names(coef(lm3FR))[-1]
  RegionFR[[r]]$lm3$model <- lm3FR
  RegionFR[[r]]$lm3$test <- TestModel(lm3FR, AfrFR)
  RegionFR[[r]]$lm3$perf <- ModelFitStats(lm3FR)
  RegionFR[[r]]$lm3$VIF <- VIF(lm3FR)
  RegionFR[[r]]$lm3$VDP <- VDP(lm3FR)
  
  # Storing 'TO' model results and metrics as 'lm3'
  cat('Storing TO model and quality check metrics as lm3\n')
  RegionTO[[r]]$lm3$terms <- names(coef(lm3TO))[-1]
  RegionTO[[r]]$lm3$model <- lm3TO
  RegionTO[[r]]$lm3$test <- TestModel(lm3TO, AfrTO)
  RegionTO[[r]]$lm3$perf <- ModelFitStats(lm3TO)
  RegionTO[[r]]$lm3$VIF <- VIF(lm3TO)
  RegionTO[[r]]$lm3$VDP <- VDP(lm3TO)
  cat('\n\n\n')
  
}

# Checking 'lm3' result and marking outliers

(r <- regions[5])

# FR...

summary(RegionFR[[r]]$lm3$model)
plot(RegionFR[[r]]$lm3$model)
RegionFR[[r]]$lm3$outliers <- integer()

# TO...

summary(RegionTO[[r]]$lm3$model)
plot(RegionTO[[r]]$lm3$model)
RegionTO[[r]]$lm3$outliers <- integer()


# NEXT: SAVE RegionFR and RegionTO to image ...
save(list = c('RegionFR','RegionTO','lookup'), file = 'ModelsRegional.RData')


