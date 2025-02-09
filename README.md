# Subtropical anticyclone impacts life-history traits of a marine top predator
Ruijiao Sun, Etienne Rouby, Henri Weimerskirch, Christophe Barbraud, Karine Delord, Caroline C. Ummenhofer, Stéphanie Jenouvrier

Please contact Ruijiao Sun (ruijiaos@ucsb.edu) for suggestions, comments, improvements, etc.

## Overview
This repository contains scripts and data to recreate the main results and figures of this study.

Please note that climate data are very large and therefore not included in this repository. However, all the monthly reanalyzed sea surface temperature, sea level pressure, and surface wind data used in this study are publicly available and can be downloaded from the ECMWF Reanalysis v5 (ERA5) https://www.ecmwf.int/en/forecasts/dataset/ecmwf-reanalysis-v5. Primary production, measured as the vertical integral of total carbon fixation between 0 and 100 m depth, can be extracted from the Japanese 55-year Reanalysis (JRA-55) https://rda.ucar.edu/datasets/d628001/.

## Folder
A brief description of the contents within each folder is provided below.

### MH_covariate
- **Mascarene_cov.m** This MATLAB code processes monthly ERA5 sea level pressure reanalysis data to derive the longitude, latitude, and strength of the Mascarene High (MH) during seasons relevant to the survival, breeding, and breeding success of wandering albatrosses.
- **MH_covariate.xlsx** This Excel file contains the Mascarene High (MH) indices utilized in this study, including the longitude and latitude of its center, as well as its strength.

### MH_Adult
This folder includes the R and JAGS code for estimating the effects of the Mascarene High on adult wandering albatrosses using Bayesian Multi-state-capture-mark-recapture model (MSCMR).
- **ENV_cov.xlsx** This dataset contains the Mascarene High (MH) covariates used in this study, including: (1) mh_lat: Latitude of the Mascarene High center; (2) mh_lon: Longitude of the Mascarene High center; (3) mh_strength: Strength of the Mascarene High.
- **202401_WA.RData** This dataset comprises the life-history records of wandering albatrosses used in this study.
  #### Survival
  This folder contains the R code for running the JAGS model to assess the effects of Mascarene High variability on the survival probability of adult wandering albatrosses. Executing the code requires the files ENV.cov.xlsx and 202401_WA.RData.
  #### Breed
  This folder contains the R code for running the JAGS model to assess the effects of Mascarene High variability on the breeding probability of adult wandering albatrosses. Executing the code requires the files ENV.cov.xlsx and 202401_WA.RData.
  #### Success
  This folder contains the R code for running the JAGS model to assess the effects of Mascarene High variability on the breeding success probability of adult wandering albatrosses. Executing the code requires the files ENV.cov.xlsx and 202401_WA.RData.
  #### Vital_rate_timeseries
  This folder contains the R code for running the JAGS model to estimate the vital rates time series of adult wandering albatrosses. Running the code requires life-history data 202401_WA.RData.
  
- **Adult_MH_vital_plot.m** This code processes the posterior_adult_ENV.mat output from the MSCMR model, which estimates the effects of the Mascarene High on vital rates, to generate Figure 1 presented in the main text.
- **ciplot.m** Function to fill color in confidence intervals.
- **invlogit.m** Function to transform data from the natural scale to the logit scale.
- **logit.m** Function to convert data to the logit scale.
- **Adult_vital_1980_2018.mat** This file contains the temporal variation in survival, breeding probability, and breeding success probability of adult wandering albatrosses from 1980 to 2018, as estimated by the MSCMR model.

  
### MH_juvenile
This folder includes the R and Nimble code for estimating the effects of the Mascarene High on juvenile wandering albatrosses using Bayesian MSCMR.
#### scripts_to_run_the_model 
This folder contains the R and Nimble code for running the MSCMR to assess the effects of Mascarene High variability on the vital rates of juvenile wandering albatrosses.

- **Juvenile_MH_vital_plot.m** This code processes the posterior_juvenile_ENV.mat output from the MSCMR model, which estimates the effects of the Mascarene High on vital rates, to generate Figure 1 presented in the main text.
- **ciplot.m** Function to fill color in confidence intervals.
- **invlogit.m** Function to transform data from the natural scale to the logit scale.
- **logit.m** Function to convert data to the logit scale.
- **Juvenile_vital_1980_2018.mat** This file contains the temporal variation in survival, breeding probability, and first-time breeding success probability of juvenile wandering albatrosses from 1980 to 2018 estimated from MSCMR models.

### Anomaly_composite_analysis
This folder contains the code for conducting anomaly composite analysis, which is used to explore the potential mechanisms through which variability in the Mascarene High influences the life-history traits of wandering albatrosses.
- **customecolormap.m** Code to generate customized color map.
- **Adult_vital_1980_2018.mat** This file contains the temporal variation in survival, breeding probability, and breeding success probability of adult wandering albatrosses from 1980 to 2018, as estimated by the JAGS model.
- **anomaly_composite_map_wind.m** This code conducts composite anomaly analysis for sea level pressure and wind conditions, generating Figure 2 as presented in the main text.
