##--------------------------------------------------------------------------------------------------------
## SCRIPT : Run wandering Albatross demographic model on real data
##
## Authors : Floriane Plard
## Last update : 2023-15-19
## R version 4.0.4 (2021-02-15) -- "Lost Library Book"
## Copyright (C) 2021 The R Foundation for Statistical Computing
## Platform: x86_64-w64-mingw32/x64 (64-bit)
##--------------------------------------------------------------------------------------------------------

# Remove all objects from the environment
rm(list = ls())

# Load the relevant packages
lapply(c("tidyverse", "dplyr", "tidyr", "nimble", "glue", "reshape2", "readxl"), library, character.only = TRUE)

# Define directory paths for data, functions, output, and model
WorkDir <- getwd() # Working directory
DataDir <- paste(WorkDir, "data", sep = "/") # Data folder
OutDir <- paste(WorkDir, "output", sep = "/") # Output folder
ModDir <- paste(WorkDir, "model", sep = "/") # Models folder
FunDir <- paste(WorkDir, "function", sep = "/") # Models folder
SourceDir <- paste(WorkDir, "source", sep = "/") # Models folder
source(glue("{FunDir}/functions.R"))

source(glue("{ModDir}/model_trends.R"))
model <- model_t
m_ <- "model_14"
sex <- "M"

# Data loading and reshape -----------------------------------------------------
# Here we keep all individuals to avoid selecting only individuals that have 
# reach a stage where they could be sexed. This would make phiPB close to 1!!!
CMR_df <- read_csv(glue("{DataDir}/data_cmr/GA_CRO_AnalyseCroisee_updated_to_2020.csv")) %>%
  mutate(across(where(is.character), ~na_if(., "NA"))) %>%
  filter(StatutBaguage == "P")

### Assign Sex and select only males
CMR_df$Sex <- NA

for (i in 1:nrow(CMR_df)) {
  # Assign the good sex given genetic 
  if (!is.na(CMR_df$SexeGenetique[i]) & CMR_df$SexeGenetique[i] == "F") {
    CMR_df$Sex[i] <- "F"
  }
  if (!is.na(CMR_df$SexeGenetique[i]) & CMR_df$SexeGenetique[i] == "M") {
    CMR_df$Sex[i] <- "M"
  }
  # Now if no genetic has been done but visual assignement has been made, assign the state
  if (is.na(CMR_df$Sex[i]) & !is.na(CMR_df$SexeVisuel[i]) & CMR_df$SexeVisuel[i] == "F") {
    CMR_df$Sex[i] <- "F"
  }
  if (is.na(CMR_df$Sex[i]) & !is.na(CMR_df$SexeVisuel[i]) & CMR_df$SexeVisuel[i] == "M") {
    CMR_df$Sex[i] <- "M"
  }
  # Now assign randomly a sex if there is no sex available. Since it will allow to better estimate survival. 
  if (is.na(CMR_df$Sex[i]) & is.na(CMR_df$SexeVisuel[i]) & is.na(CMR_df$SexeVisuel[i])) {
    s <- c("F", "M")
    CMR_df$Sex[i] <- sample(s, 1, prob = c(0.5, 0.5))
  }
}

y_tmp <- CMR_df %>%
  filter(CycleBaguage >= 1980 & CycleBaguage <= 2020) %>%
  filter(Sex == sex) %>%
  select("CycleBaguage", "1980":"2020")

y_tmp <- y_tmp[,1:42] %>%
  as.matrix()

# Covariate dataset reshape
COV.phi = rep(0, length(seq(1980,2020,1)))
COV.psi = rep(0, length(seq(1980,2020,1)))

COV.rho <- read_xlsx(paste(DataDir, "/data_cov/", "success_cov.xlsx", sep = "")) %>% 
  select(mh_m, Season, sex) %>% 
  filter(sex == sex & Season >= 1980 & Season <= 2020) %>% 
  mutate(std_mh_m = (mh_m - mean(mh_m)) / sd(mh_m)) %>% 
  pull(std_mh_m)

param_rm <- c("b.phi.PB", "b.psi.PB")

source(glue("{SourceDir}/source.R"))