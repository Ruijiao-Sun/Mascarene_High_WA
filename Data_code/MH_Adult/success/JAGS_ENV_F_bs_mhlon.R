
rm(list = ls())#clear lists
my.packs <- c("coda", "mvtnorm", "rjags", "dplyr", "reshape", "readxl",'jagsUI',"pracma")
if (any(!my.packs %in% installed.packages()[, 'Package']))install.packages(my.packs[which(!my.packs %in% installed.packages()[, 'Package'])],dependencies = TRUE,configure.args="--enable-rpath")
lapply(my.packs, require, character.only = TRUE)

load("202401_WA.RData")
cov<-read_excel("ENV_cov.xlsx",sheet="F_success")
COV<-scale(cov$mh_lon)
COV<-detrend(COV)
#only let individual enter at their first reproduction
for (i in 1:nrow(DETECTED))
{
  if (length(which(BREED[i,]==1))==0)
  {
    BREED[i,]<-NA
    SURVIVAL[i,]<-NA
    OFFSPRING[i,]<-NA
    DETECTED[i,]<-0
  }else{
         first_b <- min(which(BREED[i,]==1))
         if (AGE[i,first_b]<6){DETECTED[i,]<-0}
         else{
             BREED[i,1:(first_b-1)]<-NA
             SURVIVAL[i,1:(first_b-1)]<-NA
             OFFSPRING[i,1:(first_b-1)]<-NA
             DETECTED[i,1:(first_b-1)]<-0
             if(is.na(OFFSPRING[i,first_b])){OFFSPRING[i,first_b]<-0}
            }
        }
}
any1<-function(x){any(x==1)}
AGE<-AGE[apply(DETECTED,1,any1),]
BREED<-BREED[apply(DETECTED,1,any1),]
OFFSPRING<-OFFSPRING[apply(DETECTED,1,any1),]
SCORE<-SCORE[apply(DETECTED,1,any1),]
SEX<-SEX[apply(DETECTED,1,any1),]
SURVIVAL<-SURVIVAL[apply(DETECTED,1,any1),]
DETECTED<-DETECTED[apply(DETECTED,1,any1),]
range01 <- function(x, ...){(x - min(x, ...)) / (max(x, ...) - min(x, ...))} #define a function to scale data between 0 and 1
SCORE<-range01(SCORE,na.rm=T)

mean.score = mean(SCORE,na.rm=T)
  sd.score = sd(SCORE,na.rm=T)
 tau.score = 1/sd(SCORE,na.rm=T)^2
 cov.score = seq(min(SCORE,na.rm=T),max(SCORE,na.rm=T), length.out = 1000)
    n.pred = length(cov.score)
      pool = SCORE[-which(is.na(SCORE))]
for (i in 1:length(SCORE))
{
  if(is.na(SCORE[i]))
  {SCORE[i]<-sample(pool,1)}else{next}
}

#only select female individuals
DETECTED<-DETECTED[which(SEX==0),]
SURVIVAL<-SURVIVAL[which(SEX==0),]
BREED<-BREED[which(SEX==0),]
OFFSPRING<-OFFSPRING[which(SEX==0),]
SCORE <- SCORE[which(SEX==0)]
AGE <- AGE[which(SEX==0),]
#if only use a subset of 1000 inds data for testing
 # sub_r<-sample(1:nrow(DETECTED),100,replace=F)
 # DETECTED<-DETECTED[sub_r,]
 # SURVIVAL<-SURVIVAL[sub_r,]
 # BREED<-BREED[sub_r,]
 # OFFSPRING<-OFFSPRING[sub_r,]
 # SCORE <- SCORE[sub_r]

### for jags
### bundle data
cmr_data <- list(DETECTED = DETECTED,
                 SURVIVAL = SURVIVAL,
                 BREED = BREED,
                 OFFSPRING = OFFSPRING,
                 SCORE = SCORE,
                 AGE = AGE); 
rm(SURVIVAL, BREED, OFFSPRING, DETECTED, SCORE,AGE)
jagsdata <- vector(mode = "list", length = 0)
jagsdata$n_year <- ncol(cmr_data$DETECTED)
jagsdata$n_ind <- nrow(cmr_data$DETECTED)
jagsdata$FIRST <- apply(cmr_data$DETECTED, 1, function(row_i) { min(which(row_i == 1)) })

#jagsdata$LAST <- apply(cmr_data$SURVIVAL, 1, function(row_i) { ifelse(any(row_i == 0, na.rm = TRUE), min(which(row_i == 0)), ncol(cmr_data$SURVIVAL)) })
jagsdata$LAST <- apply(cmr_data$DETECTED, 1, function(row_i) { ifelse(max(which(row_i == 1))+5 < ncol(cmr_data$SURVIVAL), max(which(row_i == 1))+5, ncol(cmr_data$SURVIVAL)) })

jagsdata$SURVIVAL <- cmr_data$SURVIVAL
jagsdata$BREED <- cmr_data$BREED
jagsdata$OFFSPRING <- cmr_data$OFFSPRING
jagsdata$DETECTED <- cmr_data$DETECTED
jagsdata$Recruited <- c(0,rep(1,56)) #vector to distinguish post- and pre-recruiment inds
jagsdata$SCORE <- cmr_data$SCORE
jagsdata$AGE <- cmr_data$AGE

### function to initialize sampler
initfct <- function(idchain, jagsdata, perfect_detection = TRUE) {
  ### covariances
  tauA <- rgamma(4, shape = 1.5, rate = 1.5)
  A <- LA <- matrix(NA, nrow = 4, ncol = 4)
  diag(A) <- abs(rnorm(4)) * 1.5
  LA[2, 1] <- 0.25 * rnorm(1)
  LA[3, 1] <- 0.25 * rnorm(1)
  LA[4, 1] <- 0.25 * rnorm(1)
  LA[3, 2] <- 0.25 * rnorm(1)
  LA[4, 2] <- 0.25 * rnorm(1)
  LA[4, 3] <- 0.25 * rnorm(1)

  tauE <- rgamma(4, shape = 1.5, rate = 1.5)
  E <- LE <- matrix(NA, nrow = 4, ncol = 4)
  diag(E) <- abs(rnorm(4)) * 1.5
  LE[2, 1] <- 0.25 * rnorm(1)
  LE[3, 1] <- 0.25 * rnorm(1)
  LE[4, 1] <- 0.25 * rnorm(1)
  LE[3, 2] <- 0.25 * rnorm(1)
  LE[4, 2] <- 0.25 * rnorm(1)
  LE[4, 3] <- 0.25 * rnorm(1)

  init <- list(mu_phi = 1.5 * rnorm(1),
               mu_rho = 1.5 * rnorm(1),
               mu_psi = 1.5 * rnorm(1),
               mu_pi = 1.5 * rnorm(1),
               gamma = rnorm(4) * log(2) / 4,
               tauA = tauA,
               A = A,
               LA = LA,
               xi_a = replicate(4, rnorm(jagsdata$n_ind)) %*% diag(1 / sqrt(tauA)),
               tauE = tauE,
               E = E,
               LE = LE,
               xi_e = replicate(4, rnorm(jagsdata$n_year)) %*% diag(1 / sqrt(tauE))
               )
  
  if(!perfect_detection) {
    ### generate coherent values for 1 when not detected
    success <- breed <- survival <- with(jagsdata, array(NA, dim = c(n_ind, n_year)))
    # loop over individuals
    for(i in 1:jagsdata$n_ind) {
      first_det <- jagsdata$FIRST[i]
      last_det <- max(which(jagsdata$DETECTED[i, ] == 1))
      # survival first
      if(last_det + 1 < jagsdata$LAST[i]) {
        death <- sample((last_det + 1):jagsdata$LAST[i], size = 1)
      }
      else {
        death <- jagsdata$LAST[i]
      }
      survival[i, first_det:death] <- 1
      if(death != jagsdata$LAST[i]) {
        survival[i, (death + 1):jagsdata$LAST[i]] <- 0
      }
      if(first_det < jagsdata$n_year) { # check that you need to do something
        if(any(jagsdata$DETECTED[i, ] == 0, na.rm = TRUE)) {
            breed[i, (first_det):jagsdata$LAST[i]] <- rbinom(length((first_det):jagsdata$LAST[i]), size = 1, prob = 0.5)
            success[i, (first_det):jagsdata$LAST[i]] <- rbinom(length((first_det):jagsdata$LAST[i]), size = 1, prob = 0.5)
        }
      }
    }
    survival <- survival 
    breed <- breed * survival
    success <- success * breed

    init$SURVIVAL <- survival * ifelse(jagsdata$DETECTED == 0, 1, NA)
    init$BREED <- breed * ifelse(is.na(jagsdata$BREED), 1, NA)
    init$OFFSPRING <- success * ifelse(is.na(jagsdata$OFFSPRING), 1, NA)
    init$mu_p <- 1.5 * rnorm(2)
    init$unscaled_sigma2_p <- rgamma(2, 1.0, 1.0)
    init$rho_p <- rgamma(2, 1.0, 1.0)
  }
  return(init)
}


#set initial condition
IV <- initfct(jagsdata = jagsdata, perfect_detection = F)  


library(jagsUI)  # The R package which makes the interface with JAGS
 
# Bundle data
jags.data <- list(SURVIVAL = jagsdata$SURVIVAL, BREED = jagsdata$BREED, OFFSPRING = jagsdata$OFFSPRING, DETECTED = jagsdata$DETECTED, AGE = jagsdata$AGE,COV=COV,
                  FIRST = jagsdata$FIRST, LAST = jagsdata$LAST, n_ind = jagsdata$n_ind, n_year = jagsdata$n_year, 
                  Recruited = jagsdata$Recruited, max_age = max(jagsdata$AGE,na.rm = TRUE),
                  SCORE = jagsdata$SCORE, mean.score = mean.score, sd.score = sd.score, tau.score = tau.score,
                  cov.score = cov.score, n.pred = n.pred)

inits <- function(){list(SURVIVAL = IV$SURVIVAL,BREED = IV$BREED, OFFSPRING = IV$OFFSPRING)}

# Parameters monitored
parameters <- c(#female mean vital rates
                 #breeding success
                "mean.Fpi_SB","mean.Fpi_FB","mean.Fpi_PSB","mean.Fpi_PFB","mean.Fpi_NB",
                "slope",
                "epsilon"
                )
# MCMC settings
ni <- 20000      # number of iteration 
nt <- 2          # Thinning interval
nb <- 4000       # number of burnin
nc <- 3          # number of chains

# Call JAGS from R (BRT 1 min)
out1 <- jags(jags.data, inits, parameters, "WA_ENV_F_bs.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

save(out1,file = "WA_ENV_F_bs_mhlon_1980_2018.RData")











