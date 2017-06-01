#************************************************************************************************
#Project Name: PCCRC PROJECT - Bayesian Chena/Salcha Population Dynamics Model
#Creator: Curry James Cunningham, SAFS, University of Washington
#Date: 12.31.15
#
#Purpose: To estimate impact of range of environmental covariates on survival with a stage structured pop. dy. model
#
# 1) Single Population Model
# 2) Integrate bycatch removals
# 3) BMA with indicator variable
# 4) Fit to age comp. data for return
# 5) Fit to Murphy survey data
#
#
#*************************************************************************************************
#VERSIONS:
#  3: Version for meeting with Peter and Milo on Tuesday, Jan. 19, 2015. Very slow.
#  4: Converted to running jags.parallel
#  5: Collapse the beta zeros for marine stages
#  6: Single Capacity
#  7: Estimate Fmort and selectivity for bycatch
#
#ASSUMPTIONS:
#  1) Fixed proportion of offshore survey data are middle Yukon, and fixed proportions of those are Chena/Salcha.
#
#
#NOTES:
#  a) in ifelse(x,a,b) - zero is converted to false
#  b) comments count, blank lines do not. 
#  c) Likelihood functions need to reference data directly, not other variable with that value: "data ~ multi()"
#  d) Cannot output elements of objects that are never initialized within the model
#  e) The n value in dmulti(p,n) must be data and must be an integer
#  f) Dirichlet-Multinomial: Estimating the dirichWeight is necessary, errors and inability to fit abund.
#       data if pred. AC not turned into proportion and multip by dirichWeight
#     Error if pred is =0, small offsets allow to run w/o error, but not fit abund well.
#  g) Using the max(1, xxx) does not correct for issues assocaited with bycatch removals.
#  h) Use the narrow prior from King et al. on coefficients. beta.prod
#
#
#MOVING FORWARD:
#  1) Continuous catch equation for Terminal Catch (Estimate Fmorts, fit to catch data)
#  2) Simulate survival differences, by marginalizing across covariates (I've done this before.)
#  3) BSAI Bycatch Re-evaluation
#  4) BSAI Bycatch sensitivity analysis
#  5) Implement in STAN
#  6) Implement in TMB for random effects.... T
#       Time-varying natural morality, time-varying maturation. Assume hyperprior for RE's
#  7) Cross-validation (n-flold, given that there are few data)
#  8) Multipopulation model (TMB), synthetic. 
#  9) Attempting to treat Juvi as Index I=q.juvi * B 
#
#TIMINGS:
# ~ 5 min for 5e3
#
# [1] "sims: 50000 thins: 10"
# [1] "START: Fri May 20 00:25:30 2016"
# [1] "END: Fri May 20 01:35:10 2016"
#
# [1] "sims: 1e+05 thins: 50"
# [1] "START: Sun May 22 02:36:17 2016"
# [1] "END: Sun May 22 05:02:12 2016"
#
# [1] "sims: 100000 thins: 50"
# [1] "START: Mon May 23 01:06:47 2016"
# [1] "END: Mon May 23 05:26:35 2016"
#
# [1] "sims: 50000 thins: 10"
# [1] "START: Mon Jun 13 01:51:53 2016"
# [1] "END: Mon Jun 13 06:16:17 2016"
#
# [1] "sims: 10000 thins: 10"
# [1] "START: Thu Jun 16 19:51:25 2016"
# [1] "END: Thu Jun 16 20:50:23 2016"
#
# [1] "sims: 1e+05 thins: 100"
# [1] "START: Tue Nov 15 23:54:16 2016"
# [1] "END: Wed Nov 16 08:19:24 2016"

#*************************************************************************************************
require(xlsx)
require(R2jags)
require(coda)
require(ggplot2)
require(gridExtra)
require(cowplot)
require(gird)
require(mcmcplots)
require(RColorBrewer)
require(boot)
require(reshape2)
require(dclone)
require(ggmcmc)
require(BEST)
require(animation)
require(fitdistrplus)

wd <- getwd()#'/Users/curryc2/Documents/PCCRC/Jags/Final2'

# setwd(wd)

#================================================================================
##### CONTROL SECTION #####
pops <- c('chena','salcha')
n.pops <- length(pops)

#Proportion of Middle Yukon River, that is comprised of Chena/Salcha
prop.Eiler <- 0.427

#Number of Simulations and Thinning Rate:
sims <-  5e5
thins <- 500
chains <- 5

#Population to Fit:
# pop <- 'chena'
pop <- 'salcha'

run.parallel <- TRUE

#Name of Model:
# title <- 'TESTcov_SuperFinal_TotalPink_UnifPriorOE(0,5)_Gamma(1,1)_(0,1.5)PriorforBetaZero_ESTpind(2,8)_SlpPrior_N(0,0.5)'
# title <- 'FIXEDsel5_ExtremeSuperFinal_TotalPink_UnifPriorOE(0,5)_Gamma(1,1)_(0,1.5)PriorforBetaZero_ESTpind(2,8)_SlpPrior_N(0,0.5)'
# title <- 'SpawnPrior_SuperFinal_TotalPink_UnifPriorOE(0,5)_Gamma(1,1)_(0,1.5)PriorforBetaZero_ESTpind(2,8)_SlpPrior_N(0,0.5)'
# title <- 'StdSel_SuperFinal_TotalPink_UnifPriorOE(0,5)_Gamma(1,1)_(0,1.5)PriorforBetaZero_ESTpind(2,8)_SlpPrior_N(0,0.5)'
title <- 'FINAL'

model <- paste(title,' sims_',sims,' thins_',thins, sep='')

wrapper_func <- function(pop, sims=sims, thins=thins, model=model, run.parallel=TRUE) {
  
  
  #Whether or not to fit the model
  fit.model <- TRUE
  # update <- TRUE
  
  stage.abund <- c(1,4) #Spawners and after downstream/nearshore
  
  #================================================================================
  ##### READ IN DATA ######
  spawn.dat <- read.xlsx(file='Data/JAGS_ChenaSalchaData.xlsx', sheetName='Spawners')
  juvi.dat <- read.xlsx(file='Data/JAGS_ChenaSalchaData.xlsx', sheetName='Juvenile')
  juviGSI.dat <- read.xlsx(file='Data/JAGS_ChenaSalchaData.xlsx', sheetName='JuvenileGSI')
  
  catch.dat <- read.xlsx(file='Data/JAGS_ChenaSalchaData.xlsx', sheetName='Catch')
  
  bycatch.dat <- read.xlsx(file='Data/JAGS_ChenaSalchaData.xlsx', sheetName='Bycatch')
  bycatchAC.dat <- read.xlsx(file='Data/JAGS_ChenaSalchaData.xlsx', sheetName='BycatchAC', stringsAsFactors=FALSE) #Bycatch agecomp
  bycatchACsampleSize.dat <- read.xlsx(file='Data/JAGS_ChenaSalchaData.xlsx', sheetName='BycatchACsampleSize', stringsAsFactors=FALSE)
  BCgen.dat <- read.xlsx(file='Data/JAGS_ChenaSalchaData.xlsx', sheetName='BCgen')
  
  stage.dat <- read.xlsx(file='Data/JAGS_ChenaSalchaData.xlsx', sheetName='Stages', stringsAsFactors=FALSE)
  covar.dat <- read.xlsx(file='Data/JAGS_ChenaSalchaData.xlsx', sheetName='Covars', stringsAsFactors=FALSE)
  extras.dat <- read.xlsx(file='Data/JAGS_ChenaSalchaData.xlsx', sheetName='Extras')
  fecundity.dat <- read.xlsx(file='Data/JAGS_ChenaSalchaData.xlsx', sheetName='Fecundity')
  agecomp.dat <- read.xlsx(file='Data/JAGS_ChenaSalchaData.xlsx', sheetName='AgeComp')
  femAge.dat <- read.xlsx(file='Data/JAGS_ChenaSalchaData.xlsx', sheetName='FemaleAge', stringsAsFactors=FALSE)
  femProp.dat <- read.xlsx(file='Data/JAGS_ChenaSalchaData.xlsx', sheetName='PropFemale', stringsAsFactors=FALSE)
  hr.dat <- read.xlsx(file='Data/JAGS_ChenaSalchaData.xlsx', sheetName='HarvestRate', stringsAsFactors=FALSE)
  
  spawnAbundOE.dat <- read.xlsx(file='Data/JAGS_ChenaSalchaData.xlsx', sheetName='SpawnAbundOE', stringsAsFactors=FALSE)
  propEsc.dat <- read.xlsx(file='Data/JAGS_ChenaSalchaData.xlsx', sheetName='PropEsc', stringsAsFactors=FALSE) #Chena and Salcha Proportion of Escapement.
  #Extract Useful Info
  
  #Setup Population Parameters
  spawn.abund.years <- spawn.dat$year #c(1987:2014)
  juvi.abund.years <- juvi.dat$year #c(2003:2007,2009:2014)
  
  #Stages Info
  n.fw.stages <- 2
  n.o.stages <- 5
  n.stages <- n.fw.stages+n.o.stages
  stage.offset <- stage.dat$Offset #Offsets for stages
  stage.names <- stage.dat$Name #Names of stages
  
  #Brood Years to Model
  years <- spawn.abund.years - n.stages  #1980 - 2007
  n.years <- length(years)
  
  #Initial spawning abundances
  init.years <- min(years):(min(spawn.abund.years)-1)  #1980 - 1986
  n.init.years <- length(init.years)
  
  #================================================================================
  ##### COVARIATE DATA #####
  #Remove Covars that are not included
  loc.include <- which(covar.dat[2,-1]=='y')
  #Remove Upstream discharge
  if(pop=='chena') {
    loc.include <- loc.include[-2]
  }else {
    loc.include <- loc.include[-1]
  }
  
  n.covars <- length(loc.include)
  covar.names <- vector(length=n.covars)
  covar.names.full <- vector(length=n.covars)
  covar.stage <- vector(length=n.covars)
  covar.offset <- vector(length=n.covars)
  covar.link <- vector(length=n.covars)
  
  i <- 1
  for(i in 1:n.covars) {
    covar.link[i] <- covar.dat[1,loc.include[i]+1]
    covar.names[i] <- covar.dat[3,loc.include[i]+1]
    covar.names.full[i] <- covar.dat[4,loc.include[i]+1]
    covar.stage[i] <- covar.dat[5,loc.include[i]+1]
    covar.offset[i] <- as.numeric(covar.dat[6,loc.include[i]+1])
  }#next i
  
  temp.cov <- read.xlsx(file='Data/JAGS_ChenaSalchaData.xlsx', sheetName='Covars', startRow=7)
  covar.years <- temp.cov[,1]
  covars <- temp.cov[,loc.include+1]
  
  #Determine number of prod and cap covariates
  loc.prod.covars <- which(covar.link=='prod')
  n.prod.covars <- length(loc.prod.covars)
  prod.covar.names <- covar.names[loc.prod.covars]
  prod.covar.names.full <- covar.names.full[loc.prod.covars]
  
  loc.cap.covars <- which(covar.link=='cap')
  n.cap.covars <- length(loc.cap.covars)
  cap.covar.names <- covar.names[loc.cap.covars]
  cap.covar.names.full <- covar.names.full[loc.cap.covars]
  
  #================================================================================
  ##### ABUNDANCE DATA #####
  
  #ADULT ABUNDANCE
  spawn.abund <- spawn.dat[,which(names(spawn.dat)==pop)]
  
  #JUVENILE ABUNDANCE
  # juvi.abund <- juvi.dat[,which(names(juvi.dat)==pop)]
  juvi.idx <- juvi.dat$index
  n.juvi.abund <- length(juvi.abund.years)
  juvi.abund <- length(n.juvi.abund) #Stock-specific juvenile abundance
  
  y <- 1
  for(y in 1:n.juvi.abund) {
    year <- juvi.abund.years[y]
    loc.yr <- which(juviGSI.dat$year==year)
    
    #Calculate MYR idx based on GSI data from Murphy
    if(length(loc.yr)==0) {  #Genetic data are NOT available
      juvi.abund[y] <- juvi.idx[y] * mean(juviGSI.dat$mean)
    }else {  #Genetic data are available
      juvi.abund[y] <- juvi.idx[y] * juviGSI.dat$mean[loc.yr]
    }
  }#next y
  
  #Translade from MYR to Chena/Salcha - Based on Eiler tracking data
  juvi.abund <- juvi.abund * prop.Eiler
  
  #Translate to index for the specific stock
  juvi.abund <- juvi.abund * mean(propEsc.dat$prop[propEsc.dat$pop==pop])
  
  #Establish years for likelihood in jags code
  temp.years <- (juvi.abund.years-stage.offset[2]) #Brood years
  
  #Brood years fit
  juvi.broodYears.fit <- temp.years[which(temp.years <= max(years))]
  
  #Location in Brood years for N matrix
  loc.juvi.N <-  which(years %in% juvi.broodYears.fit) #Location in brood years
  
  #Calendar years fit
  juvi.calendarYears.fit <- juvi.broodYears.fit + stage.offset[2]
  
  #Location in juvenile index
  loc.juvi.abund <- which(juvi.abund.years %in% juvi.calendarYears.fit)
  
  #Number of likelihood comparisons
  n.juvi.locs <- length(loc.juvi.N)
  #================================================================================
  ##### FECUNDITY DATA #####
  #Deterine fecundity for each age class (2-5)
  fecundity <- fecundity.dat[,which(names(fecundity.dat)==pop)]
  fec.age <- fecundity.dat$OceanAge
  
  #================================================================================
  ##### AGE COMPOSITION DATA #####
  #Establish Correct Agecomp Data
  # agecomp <- as.matrix(agecomp.dat[agecomp.dat$stock==pop,-c(1:2)]) #Fudge Factor
  
  agecomp <-  round(agecomp.dat[agecomp.dat$stock==pop,-c(1:2)],0)
  oAges <- 1:5
  n.oAges <- length(oAges)
  total <- apply(agecomp,1,sum) #Total for the multinom likelihood
  
  #Age composition years
  agecomp.years <- agecomp.dat$year[agecomp.dat$stock==pop]
  n.agecomp.years <- length(agecomp.years)
  
  #Female age
  fem.oAges <- 2:5
  n.fem.oAges <- length(fem.oAges)
  
  #Female Age Composition Data
  femAge <- femAge.dat[femAge.dat$stock==pop,-c(1:2)]
  
  #Female Spawning Proportion
  femProp <- femProp.dat[femProp.dat$stock==pop, -c(1:2)]
  
  #================================================================================
  ##### COVARIATE REFERENCES #####
  
  prod.ref <- matrix(data=0, nrow=n.stages, ncol=n.covars)
  n.prod.ref <- vector(length=n.stages)
  
  cap.ref <- matrix(data=0, nrow=n.stages, ncol=n.covars)
  n.cap.ref <- vector(length=n.stages)
  
  s <- 1
  for(s in 1:n.stages) {
    temp.stage <- stage.names[s]
    #Prod
    temp.loc <- which(covar.stage==temp.stage & covar.link=='prod')
    n.temp.loc <- length(temp.loc)
    n.prod.ref[s] <- n.temp.loc
    if(n.temp.loc > 0) {
      prod.ref[s,temp.loc] <- temp.loc
    }
    #Cap
    temp.loc <- which(covar.stage==temp.stage & covar.link=='cap')
    n.temp.loc <- length(temp.loc)
    n.cap.ref[s] <- n.temp.loc
    if(n.temp.loc > 0) {
      cap.ref[s,temp.loc] <- temp.loc
    }
  }#next s
  
  #================================================================================
  ##### BYCATCH DATA ##### 
  #Calculate bycatch
  bycatch.years <- bycatch.dat$year
  n.bycatch.years <- length(bycatch.years)
  
  #Annual bycatch data
  bycatch <- vector(length=n.bycatch.years)
  
  #Calculate MYR Bycatch
  y <- 1
  for(y in 1:n.bycatch.years) {
    temp.year <- bycatch.years[y]
    
    temp.loc <- which(BCgen.dat$year==temp.year)
    if(length(temp.loc)==0) {
      bycatch[y] <- bycatch.dat$bycatch[y] * mean(BCgen.dat$prop)
    }else {
      bycatch[y] <- bycatch.dat$bycatch[y] * BCgen.dat$prop[temp.loc]
    }
  }#next y
  
  #Translade from MYR to Chena/Salcha - Based on Eiler tracking data
  bycatch <- bycatch * prop.Eiler
  
  #Translate to index for the specific stock
  bycatch <- bycatch * mean(propEsc.dat$prop[propEsc.dat$pop==pop])
  
  #================================================================================
  ##### BYCATCH AGECOMP DATA ##### - continue here!
  #Bycatch Age Composition
  ac.bycatch.years <- bycatchAC.dat$year
  n.ac.bycatch.years <- length(ac.bycatch.years)
  #Data object for bycatch agecomp
  ac.bycatch <- bycatchAC.dat[,-1]
  
  y <- 1
  for(y in 1:n.ac.bycatch.years) {
    ac.bycatch[y,] <- ac.bycatch[y,]/sum(ac.bycatch[y,])
    ac.bycatch[y,] <- round(ac.bycatch[y,]*(bycatchACsampleSize.dat$sampleSize[y]),0)
  }#next y
  
  #Multiply by sample size
  sum_ac.bycatch <- apply(ac.bycatch,1,sum)#bycatchACsampleSize$sampleSize
  
  #Fishing mortality
  n.fmort <- n.years+max(stage.offset)-stage.offset[n.fw.stages+1] #Fmort = 1983 - 2014
  
  #Catch Data
  catch.years <- catch.dat$year
  catch <- catch.dat[,which(names(catch.dat)==pop)]
  
  #Harvest Rate Data
  hr <- hr.dat[,which(names(hr.dat)==pop)]
  
  #Determine Prior for Spawning Abundance OE
  spawnAbundOE <- spawnAbundOE.dat[,which(names(spawnAbundOE.dat)==pop)]
  #Remove NA's
  spawnAbundOE <- spawnAbundOE[-which(is.na(spawnAbundOE))]
  #Fit lognormal distribution to get prior parameters
  fit.oe <- fitdist(spawnAbundOE, distr='lnorm', method='mle')
  #Prior specifications
  meanlog.spawnOE <- as.numeric(fit.oe$estimate[1])
  sdlog.spawnOE <- as.numeric(fit.oe$estimate[2])
  #Determine necessary covariate years
  # min.offset <- 0
  # max.offset <- 3
  # covar.years <- (min(years)+min.offset):(max(years)+max.offset)
  
  #================================================================================================================
  yukon_bayes <- NULL
  ##### JAGS MODEL CODE ######
  yukon_bayes <- function() {
    
    ### PRIORS ###
    # for(s in 1:n.stages) {
    for(s in 1:(n.fw.stages+1)) {
     
      beta.zero.prod[s] ~ dnorm(0,pow(1.5,-2)) 
      # beta.zero.prod[s] ~ dnorm(0,pow(5,-2)) 
#       exp.beta.zero.cap[s] ~ dunif(1e3,5e7)
#       beta.zero.cap[s] <- log(exp.beta.zero.cap[s])
    
    }#next s
    #Add in capacity for additional stages
    slp.beta.zero.prod ~ dnorm(0,pow(0.5,-2))
    #Gelman (2008) Recommendation: Cauchy, center 0, scale 2.5, Note: Cauchy is special case of t-dist with df=1
    # slp.beta.zero.prod ~ dt(0,2.5,1)
    
    #Change for each successive year
    for(s in 2:n.o.stages) {
      beta.zero.prod[s+n.fw.stages] <- beta.zero.prod[n.fw.stages+1] + (s-1)*slp.beta.zero.prod
    }
#     #Est FW fix OC capacity
#     exp.beta.zero.cap[1] ~ dunif(1e3,1e7)
#     exp.beta.zero.cap[2] <- 1e7
    
    for(s in 1:2) {
    
      #Separate Freshwater and Marine Capacities
      exp.beta.zero.cap[s] ~ dunif(1e3,1e7)
      beta.zero.cap[s] <- log(exp.beta.zero.cap[s])
    }#next s
    

   
    #Indicator variable probability
    # pind ~ dbeta(1,1)
    pind ~ dbeta(2,8)
    # pind <- 0.25
    
    beta.prior.tau <- 1/beta.prior.var #pow(beta.prior.sd,-2)

    #Heirarchical Prior for Lindley's Paradox
    sigma.prod  ~ dgamma(1,1) #dunif(0,10)
    tau.prod <- pow(sigma.prod,-2)
    
    for(p in 1:n.prod.covars) {
      # beta.prod[p] ~ dnorm(0, beta.prior.tau)
      beta.prod[p] ~ dnorm(0, tau.prod)
      ind.prod[p] ~ dbern(pind)
      # ind.prod[p] ~ dbern(0.3)
    }
    
#     for(c in 1:n.cap.covars) {
#       beta.cap[c] ~ dnorm(0, beta.prior.tau)
#       # ind.cap[c] ~ dbern(pind)
#       ind.cap[c] ~ dbern(0.3)
#     }
#     
    # init.spawners <- rep(log(mean(spawn.abund)), n.init.years)
    
    # y <- 1
    #for(y in 1:n.init.years) {
    #  init.spawners[y] ~ dnorm(6035, pow(3302,-2))
    #}#next y
    
    # mat.param <- 3 #Parameter for maturation dunif(0,7)
    mat.param ~ dunif(0,n.oAges)
    
    #Informative Prior
    # sigma.spawn ~ dlnorm(meanlog.spawnOE,pow(sdlog.spawnOE,-2)) #2
    # tau.spawn <- pow(sigma.spawn, -2)
    
    #Gamma Prior
    # tau.spawn ~ dgamma(1,0.01)#dgamma(0.001,0.001)
    # sigma.spawn <- 1/sqrt(tau.spawn)
    
    sigma.spawn ~ dunif(0,5)#dgamma(1,1)
    tau.spawn <- pow(sigma.spawn,-2)

    # sigma.juvi ~ dunif(1e-3,5e5)
#     sigma.juvi ~ dunif(1e-3,5) #2
#     tau.juvi <- pow(sigma.juvi, -2)
#     tau.juvi ~ dgamma(1,0.01)#dgamma(0.001,0.001)
#     sigma.juvi <- 1/sqrt(tau.juvi)
    
    sigma.juvi ~ dunif(0,5)#dgamma(1,1)
    tau.juvi <- pow(sigma.juvi,-2)
    
    #Gamma Prior
#     tau.spawn ~ dgamma(0.1,0.1)
#     sigma.spawn <- sqrt(1/tau.spawn)
#     
#     tau.juvi ~ dgamma(0.1,0.1)
#     sigma.juvi <- sqrt(1/tau.juvi)
    
    #FISHING MORTALITY RATE FOR BYCATCH EFFECT
    #1:32 checks out
#     for(y in (1+stage.offset[1+n.fw.stages]):(n.years+stage.offset[n.o.stages+n.fw.stages])) {
#       F_bycatch[y-stage.offset[1+n.fw.stages]] ~ dunif(0,5) #dlnorm(log(mu.Fmort), pow(sd.Fmort,-2))
#       # F_bycatch[y-stage.offset[1+n.fw.stages]] ~ dunif(0,1)
#     }#next y
    
    #HEIRARCHICAL BYCATCH FISHING MORTALITY RATE
    mu.Fmort ~ dunif(0,4)
    sd.Fmort ~ dunif(1e-3,2)
    for(y in (1+stage.offset[1+n.fw.stages]):(n.years+stage.offset[n.o.stages+n.fw.stages])) {
      F_bycatch[y-stage.offset[1+n.fw.stages]] ~ dlnorm(log(mu.Fmort), pow(sd.Fmort,-2))
      #Translate into Harvest Rate
      # HR_bycatch[y-stage.offset[1+n.fw.stages]] <- 1-exp(-1*F_bycatch[y-stage.offset[1+n.fw.stages]])
    }#next y

    #INDEPENDENT SELECTIVITY ORIGINAL
    # for(a in 1:n.o.stages) {
    #  temp_sel_bycatch[a] ~ dunif(0,1)
    # }#next a
    # 
    # for(a in 1:n.o.stages) {
    #   sel_bycatch[a] <- temp_sel_bycatch[a]/max(temp_sel_bycatch)
    # }
    #FIX 5 SELECTIVITY equal to 4
    for(a in 1:(n.o.stages-1)) {
      sel_bycatch[a] ~ dunif(0,1)
    }#next a
    sel_bycatch[n.o.stages] <- sel_bycatch[n.o.stages-1]
    
    #LOGISTIC SELECTIVITY
    
    # selex_1 ~ dunif(0.001,20)#1, n.o.stages)  #Age at inflection
    # selex_2 ~ dunif(0.001,20)  #Width for 95% selection
    # for(a in 1:n.o.stages) {
    #   sel_bycatch[a] <- 1/(1+exp(-1*log(19) *(a-selex_1)/selex_2 ))
    # }#next a
    
    #Single 1-ocean selectivity
#     sel_bycatch[1] ~ dunif(0,1)
#     sel_bycatch[2] ~ dunif(0,1)
#     sel_bycatch[3] ~ dunif(0,1)
#     sel_bycatch[4] ~ dunif(0,1)
#     sel_bycatch[5] <- 1
    
#     sigma.bycatch ~ dunif(1e-3,5)
#     tau.bycatch <- pow(sigma.bycatch, -2)
    #Gamma
#     tau.bycatch ~ dgamma(1,0.01)#dgamma(0.001,0.001)
#     sigma.bycatch <- 1/sqrt(tau.bycatch)
    
    sigma.bycatch ~ dunif(0,5)#dgamma(1,1)
    tau.bycatch <- pow(sigma.bycatch,-2)
    
#     tau.bycatch ~ dgamma(1,0.001)
#     sigma.bycatch <- sqrt(1/tau.bycatch)
    
    ### CREATE DATA OBJECTS ###
    ### CREATE MATURATION SCHEDULE ###
    # a <- 1
    for(a in 1:n.oAges) {
      prob.mat[a] <- 1/(1+exp(mat.param*(5-a)+log(1/0.99-1))) # Probability of maturing at age
    }
    # prob.mat[n.oAges] <- 1
    
    ### JUVENILE CATCHABILITY ###
    # q.juvi ~ dunif(0, 10) #I=qB
    
    ### CREATE COVARIATE MATRIX ###
    
    #   prod.mtx[1,prod.ref[1,1:n.prod.ref[1]]] <- prod.ref[1,1:n.prod.ref[1]]
    #   prod.mtx[1,(prod.ref[n.prod.ref[1]]+1):n.covars] <- prod.ref[1,prod.ref[1,(prod.ref[n.prod.ref[1]]+1):n.covars]]
    
    #PRODUCTIVITY
    #Stage 1
    prod.mtx[1,1:4] <- ind.prod[prod.ref[1,1:4]]*beta.prod[prod.ref[1,1:4]] #Covars
    prod.mtx[1,5:n.covars] <- prod.ref[1,5:n.covars] 
    #Stage 2
    prod.mtx[2,1:4] <- prod.ref[2,1:4]
    prod.mtx[2,5:10] <- ind.prod[prod.ref[2,5:10]]*beta.prod[prod.ref[2,5:10]] #Covars
    prod.mtx[2,11:n.covars] <- prod.ref[2,11:n.covars]
    #Stage 3
    prod.mtx[3,1:10] <- prod.ref[3,1:10]
    prod.mtx[3,11:n.covars] <- ind.prod[prod.ref[3,11:n.covars]]*beta.prod[prod.ref[3,11:n.covars]] #Covars
    # prod.mtx[3,13:n.covars] <- prod.ref[3,13:n.covars]
    #Stage4-7
    prod.mtx[4,1:n.covars] <- prod.ref[4,1:n.covars]
    prod.mtx[5,1:n.covars] <- prod.ref[5,1:n.covars]
    prod.mtx[6,1:n.covars] <- prod.ref[6,1:n.covars]
    prod.mtx[7,1:n.covars] <- prod.ref[7,1:n.covars]
    
    #CAPACITY
    #Stage 1
#     cap.mtx[1,1:n.covars] <- cap.ref[1,1:n.covars]
#     #Stage 2
#     cap.mtx[2,1:n.covars] <- cap.ref[2,1:n.covars]
#     #Stage 3
#     cap.mtx[3,1:12] <- cap.ref[3,1:12] 
#     cap.mtx[3,13:n.covars] <- ind.cap[1:n.cap.covars]*beta.cap[1:n.cap.covars] #Covars
    #Stage 4-7
    cap.mtx[1,1:n.covars] <- cap.ref[1,1:n.covars]
    cap.mtx[2,1:n.covars] <- cap.ref[2,1:n.covars]
    cap.mtx[3,1:n.covars] <- cap.ref[3,1:n.covars]
    cap.mtx[4,1:n.covars] <- cap.ref[4,1:n.covars]
    cap.mtx[5,1:n.covars] <- cap.ref[5,1:n.covars]
    cap.mtx[6,1:n.covars] <- cap.ref[6,1:n.covars]
    cap.mtx[7,1:n.covars] <- cap.ref[7,1:n.covars]
    
    ### SIMULATE LIFECYCLE ###
    ### Initial Years ###
    for(y in 1:n.init.years) {
      #Spawning Abundance
      est.init.spawners[y] ~ dunif(1,1e4)#dnorm(mean(spawn.abund),pow(sd(spawn.abund),-2));I(1,1e6)#(1,1e4)
      N[y,1] <- est.init.spawners[y]  #ifelse(y<=3, init.spawners[y], 5)
      
      #Eggs
      for(o in 1:n.fem.oAges) {
        EpF[y,o] <- femAge[y,o]*fecundity[o] #Eggs-per-female
      }#next o
      #Eggs[y] <- femProp[y]*sum(EpF[y,])
      N[y,2] <- N[y,1]*femProp[y]*sum(EpF[y,])
      
      #Loop Freshwater Stages
      for(s in 1:n.fw.stages) {
        #Calculate Covariate Effects
        for(c in 1:n.covars) {
          prod.eff[y,s,c] <- prod.mtx[s,c] * covars[y+covar.offset[c],c] #Works
          cap.eff[y,s,c] <- cap.mtx[s,c] * covars[y+covar.offset[c],c]
        }
        
        #Calculate Stage Prod/Cap
        #sum.prod.eff[y,1] <- sum(prod.eff[y,1,])#Gives same answer=GOOD!
        prod[y,s] <- 1/(1+exp(-1*beta.zero.prod[s] - sum(prod.eff[y,s,])))
        cap[y,s] <- exp(beta.zero.cap[1] + sum(cap.eff[y,s,]))
        
        #Calculate Survival
        surv[y,s] <- prod[y,s]/(1+((prod[y,s]*N[y,s+1])/cap[y,s])) 
        #Update Pop Matrix
        N[y,s+2] <- N[y,s+1] * surv[y,s]
        
      }#next fw stage (s)
      
      #Loop Ocean Stages
      for(s in 1:n.o.stages) {
        #Calculate Covariate Effects
        for(c in 1:n.covars) {
          prod.eff[y,s+n.fw.stages,c] <- prod.mtx[s+n.fw.stages,c] * covars[y+covar.offset[c],c] #Works
          cap.eff[y,s+n.fw.stages,c] <- cap.mtx[s+n.fw.stages,c] * covars[y+covar.offset[c],c]
        }
        
        #Calculate Stage Prod/Cap
        #sum.prod.eff[y,1] <- sum(prod.eff[y,1,])#Gives same answer=GOOD!
        #         prod[y,s+n.fw.stages] <- 1/(1+exp(-1*beta.zero.prod[s+n.fw.stages] - sum(prod.eff[y,s+n.fw.stages,])))
        #         cap[y,s+n.fw.stages] <- exp(beta.zero.cap[s+n.fw.stages] + sum(cap.eff[y,s+n.fw.stages,]))
        
        prod[y,s+n.fw.stages] <- 1/(1+exp(-1*beta.zero.prod[s+n.fw.stages] - sum(prod.eff[y,s+n.fw.stages,])))
        cap[y,s+n.fw.stages] <- exp(beta.zero.cap[2] + sum(cap.eff[y,s+n.fw.stages,]))
        
        #Calculate Survival
        surv[y,s+n.fw.stages] <- prod[y,s+n.fw.stages]/(1+((prod[y,s+n.fw.stages]*N[y,(s+n.fw.stages+1)])/cap[y,s+n.fw.stages])) 
        
        #Calculate bycatch
#         bycatch.mtx[y+stage.offset[s+n.fw.stages]-stage.offset[1+n.fw.stages],s] <- N[y,s+n.fw.stages+1] * surv[y,s+n.fw.stages]*
#                                   (1-exp(-1*sel_bycatch[s]*F_bycatch[y+stage.offset[s+n.fw.stages]-stage.offset[1+n.fw.stages]]))
        
        #Continuous
        bycatch.mtx[y+stage.offset[s+n.fw.stages]-stage.offset[1+n.fw.stages],s] <- ((sel_bycatch[s]*F_bycatch[y+stage.offset[s+n.fw.stages]-stage.offset[1+n.fw.stages]])/
                                                                                      (sel_bycatch[s]*F_bycatch[y+stage.offset[s+n.fw.stages]-stage.offset[1+n.fw.stages]] + 
                                                                                         -1*log(surv[y,s+n.fw.stages]))) *
                                                                                    N[y,s+n.fw.stages+1] *
                                                                                    (1-exp(-1*( (sel_bycatch[s]*F_bycatch[y+stage.offset[s+n.fw.stages]-stage.offset[1+n.fw.stages]])+ 
                                                                                                  (-1*log(surv[y,s+n.fw.stages])) )))
        
        #As Harvest Rate
#         bycatch.mtx[y+stage.offset[s+n.fw.stages]-stage.offset[1+n.fw.stages],s] <- N[y,s+n.fw.stages+1] * surv[y,s+n.fw.stages]*
#                                    ((F_bycatch[y+stage.offset[s+n.fw.stages]-stage.offset[1+n.fw.stages]]))
        
        #Update Pop Matrix
#         N[y,s+n.fw.stages+2] <- (N[y,s+n.fw.stages+1] * surv[y,s+n.fw.stages] - bycatch.mtx[y+stage.offset[s+n.fw.stages]-stage.offset[1+n.fw.stages],s]) * 
#                                   (1-prob.mat[s+n.fw.stages-2])
        
        #Continuous
        N[y,s+n.fw.stages+2] <- (N[y,s+n.fw.stages+1] * exp(-1*( (-1*log(surv[y,s+n.fw.stages])) + (sel_bycatch[s]*F_bycatch[y+stage.offset[s+n.fw.stages]-stage.offset[1+n.fw.stages]]) ))) * 
          (1-prob.mat[s+n.fw.stages-2])
        
        #Bycatch First
#         N[y,s+n.fw.stages+2] <- (N[y,s+n.fw.stages+1]  - bycatch.mtx[y+stage.offset[s+n.fw.stages]-stage.offset[1+n.fw.stages],s]) * surv[y,s+n.fw.stages] * 
#           (1-prob.mat[s+n.fw.stages-2])
        
        

        #Returns
#         Ret[y,s] <- (N[y,s+n.fw.stages+1] * surv[y,s+n.fw.stages] - bycatch.mtx[y+stage.offset[s+n.fw.stages]-stage.offset[1+n.fw.stages],s]) * 
#                       (prob.mat[s+n.fw.stages-2])
        
        #Continuous
        Ret[y,s] <- (N[y,s+n.fw.stages+1] * exp(-1*( (-1*log(surv[y,s+n.fw.stages])) + (sel_bycatch[s]*F_bycatch[y+stage.offset[s+n.fw.stages]-stage.offset[1+n.fw.stages]]) ))) * 
          (prob.mat[s+n.fw.stages-2])
        
        #Bycatch First
#         Ret[y,s] <- (N[y,s+n.fw.stages+1]  - bycatch.mtx[y+stage.offset[s+n.fw.stages]-stage.offset[1+n.fw.stages],s]) * surv[y,s+n.fw.stages] * 
#           (prob.mat[s+n.fw.stages-2])
        
        # Ret[y,s] <- (N[y,s+n.fw.stages+1] * surv[y,s+n.fw.stages] - bycatch[y+stage.offset[s+n.fw.stages],s+n.fw.stages-2])*(prob.mat[s+n.fw.stages-2])
        
        #     Catch[y,s] <- Ret[y,s] * (1-exp(-1*fmort[y+stage.offset[s]-3]))
        
        #     Esc[y,s] <- Ret[y,s] * exp(-1*fmort[y+stage.offset[s]-3])
        
        Catch[y,s] <- Ret[y,s]*hr[y+stage.offset[s+n.fw.stages]]
        
        Esc[y,s] <- Ret[y,s]*(1-hr[y+stage.offset[s+n.fw.stages]])
        
        Esc.mtx[y+stage.offset[s+n.fw.stages],s] <- Esc[y,s]
        
        # N[y+stage.offset[s+n.fw.stages],1] <- N[y+stage.offset[s+n.fw.stages],1] + Esc[y,s+n.fw.stages]
      }#next o stage (s)
    }#next init.year
    
    #Fill in blanks in Esc.mtx which are not filled
    #   Esc.mtx[1,2] <- 0
    #   Esc.mtx[1,3] <- 0
    #   Esc.mtx[1,4] <- 0
    #   Esc.mtx[1,5] <- 0
    #   
    #   Esc.mtx[2,3] <- 0
    #   Esc.mtx[2,4] <- 0
    #   Esc.mtx[2,5] <- 0
    #   
    #   Esc.mtx[3,4] <- 0
    #   Esc.mtx[3,5] <- 0
    #   
    #   Esc.mtx[4,5] <- 0
    
    ### Main Years ###
    for(y in (1+n.init.years):n.years) { #double check this
      #Spawning Abundance
      
      N[y,1] <- sum(Esc.mtx[y,])  #ifelse(y<=3, init.spawners[y], 5)
      
      #Eggs
      for(o in 1:n.fem.oAges) {
        EpF[y,o] <- femAge[y,o]*fecundity[o] #Eggs-per-female
      }#next o
      #Eggs[y] <- femProp[y]*sum(EpF[y,])
      N[y,2] <- N[y,1]*femProp[y]*sum(EpF[y,])
      
      #Loop Freshwater Stages
      for(s in 1:n.fw.stages) {
        #Calculate Covariate Effects
        for(c in 1:n.covars) {
          prod.eff[y,s,c] <- prod.mtx[s,c] * covars[y+covar.offset[c],c] #Works
          cap.eff[y,s,c] <- cap.mtx[s,c] * covars[y+covar.offset[c],c]
        }
        
        #Calculate Stage Prod/Cap
        #sum.prod.eff[y,1] <- sum(prod.eff[y,1,])#Gives same answer=GOOD!
        prod[y,s] <- 1/(1+exp(-1*beta.zero.prod[s] - sum(prod.eff[y,s,])))
        cap[y,s] <- exp(beta.zero.cap[1] + sum(cap.eff[y,s,]))
        
        #Calculate Survival
        surv[y,s] <- prod[y,s]/(1+((prod[y,s]*N[y,s+1])/cap[y,s])) 
        #Update Pop Matrix
        N[y,s+2] <- N[y,s+1] * surv[y,s]
        
      }#next fw stage (s)
      
      #Loop Ocean Stages
      for(s in 1:n.o.stages) {
        #Calculate Covariate Effects
        for(c in 1:n.covars) {
          prod.eff[y,s+n.fw.stages,c] <- prod.mtx[s+n.fw.stages,c] * covars[y+covar.offset[c],c] #Works
          cap.eff[y,s+n.fw.stages,c] <- cap.mtx[s+n.fw.stages,c] * covars[y+covar.offset[c],c]
        }
        
        #Calculate Stage Prod/Cap
        #sum.prod.eff[y,1] <- sum(prod.eff[y,1,])#Gives same answer=GOOD!
        #         prod[y,s+n.fw.stages] <- 1/(1+exp(-1*beta.zero.prod[s+n.fw.stages] - sum(prod.eff[y,s+n.fw.stages,])))
        #         cap[y,s+n.fw.stages] <- exp(beta.zero.cap[s+n.fw.stages] + sum(cap.eff[y,s+n.fw.stages,]))
        
        prod[y,s+n.fw.stages] <- 1/(1+exp(-1*beta.zero.prod[s+n.fw.stages] - sum(prod.eff[y,s+n.fw.stages,])))
        cap[y,s+n.fw.stages] <- exp(beta.zero.cap[2] + sum(cap.eff[y,s+n.fw.stages,]))
        
        #Calculate Survival
        
        surv[y,s+n.fw.stages] <- prod[y,s+n.fw.stages]/(1+((prod[y,s+n.fw.stages]*N[y,(s+n.fw.stages+1)])/cap[y,s+n.fw.stages])) 
        
        # surv[y,s+n.fw.stages] <- prod[y,s+n.fw.stages]/(1+((prod[y,s+n.fw.stages]*N[y,4])/cap[y,s+n.fw.stages])) 
        
        #Calculate Bycatch
#         bycatch.mtx[y+stage.offset[s+n.fw.stages]-stage.offset[1+n.fw.stages],s] <- N[y,s+n.fw.stages+1] * surv[y,s+n.fw.stages] *
#                           (1-exp(-1*sel_bycatch[s]*F_bycatch[y+stage.offset[s+n.fw.stages]-stage.offset[1+n.fw.stages]]))

        #Continuous
        bycatch.mtx[y+stage.offset[s+n.fw.stages]-stage.offset[1+n.fw.stages],s] <- ((sel_bycatch[s]*F_bycatch[y+stage.offset[s+n.fw.stages]-stage.offset[1+n.fw.stages]])/
                                                                                       (sel_bycatch[s]*F_bycatch[y+stage.offset[s+n.fw.stages]-stage.offset[1+n.fw.stages]] + 
                                                                                          -1*log(surv[y,s+n.fw.stages]))) *
          N[y,s+n.fw.stages+1] *
          (1-exp(-1*( (sel_bycatch[s]*F_bycatch[y+stage.offset[s+n.fw.stages]-stage.offset[1+n.fw.stages]])+ 
                        (-1*log(surv[y,s+n.fw.stages])) )))
        
        #As Harvest Rate
        #         bycatch.mtx[y+stage.offset[s+n.fw.stages]-stage.offset[1+n.fw.stages],s] <- N[y,s+n.fw.stages+1] * surv[y,s+n.fw.stages]*
        #                                    ((F_bycatch[y+stage.offset[s+n.fw.stages]-stage.offset[1+n.fw.stages]]))
        
        #Update Pop Matrix
        #         N[y,s+n.fw.stages+2] <- (N[y,s+n.fw.stages+1] * surv[y,s+n.fw.stages] - bycatch.mtx[y+stage.offset[s+n.fw.stages]-stage.offset[1+n.fw.stages],s]) * 
        #                                   (1-prob.mat[s+n.fw.stages-2])
        
        #Continuous
        N[y,s+n.fw.stages+2] <- (N[y,s+n.fw.stages+1] * exp(-1*( (-1*log(surv[y,s+n.fw.stages])) + (sel_bycatch[s]*F_bycatch[y+stage.offset[s+n.fw.stages]-stage.offset[1+n.fw.stages]]) ))) * 
          (1-prob.mat[s+n.fw.stages-2])
        
        #Bycatch First
        #         N[y,s+n.fw.stages+2] <- (N[y,s+n.fw.stages+1]  - bycatch.mtx[y+stage.offset[s+n.fw.stages]-stage.offset[1+n.fw.stages],s]) * surv[y,s+n.fw.stages] * 
        #           (1-prob.mat[s+n.fw.stages-2])
        
        
        
        #Returns
        #         Ret[y,s] <- (N[y,s+n.fw.stages+1] * surv[y,s+n.fw.stages] - bycatch.mtx[y+stage.offset[s+n.fw.stages]-stage.offset[1+n.fw.stages],s]) * 
        #                       (prob.mat[s+n.fw.stages-2])
        
        #Continuous
        Ret[y,s] <- (N[y,s+n.fw.stages+1] * exp(-1*( (-1*log(surv[y,s+n.fw.stages])) + (sel_bycatch[s]*F_bycatch[y+stage.offset[s+n.fw.stages]-stage.offset[1+n.fw.stages]]) ))) * 
          (prob.mat[s+n.fw.stages-2])
        
        #Bycatch First
#         Ret[y,s] <- (N[y,s+n.fw.stages+1]  - bycatch.mtx[y+stage.offset[s+n.fw.stages]-stage.offset[1+n.fw.stages],s]) * surv[y,s+n.fw.stages] *
#           (prob.mat[s+n.fw.stages-2])
        
        # Ret[y,s] <- (N[y,s+n.fw.stages+1] * surv[y,s+n.fw.stages] - bycatch[y+stage.offset[s+n.fw.stages],s+n.fw.stages-2])*(prob.mat[s+n.fw.stages-2])
        
        #     Catch[y,s] <- Ret[y,s] * (1-exp(-1*fmort[y+stage.offset[s]-3]))
        
        #     Esc[y,s] <- Ret[y,s] * exp(-1*fmort[y+stage.offset[s]-3])
        
        Catch[y,s] <- Ret[y,s]*hr[y+stage.offset[s+n.fw.stages]]
        
        Esc[y,s] <- Ret[y,s]*(1-hr[y+stage.offset[s+n.fw.stages]])
        
        Esc.mtx[y+stage.offset[s+n.fw.stages],s] <- Esc[y,s] #Check for error
        
        # N[y+stage.offset[s+n.fw.stages],1] <- N[y+stage.offset[s+n.fw.stages],1] + Esc[y,s+n.fw.stages]
      }#next o stage (s)
    }#next init.year
    
    #=======================================================================
    #LIKELIHOOD FOR SPAWNING ABUNDANCE - good 
    for(i in (n.init.years+1):(n.years+3)) {
      #Extract Predictions
      spawn.pred[i-n.init.years] <- sum(Esc.mtx[i,])
      #Posterior predictive
      ln.spawn.post.pred[i-n.init.years] ~ dnorm(log(spawn.pred[i-n.init.years]), tau.spawn)
      spawn.post.pred[i-n.init.years] <- exp(ln.spawn.post.pred[i-n.init.years])
      #Extract Observations
      # spawn.obs[i-n.init.years] <- spawn.abund[i-n.init.years] #For plotting purpose
      spawn.obs[i-n.init.years] <- exp(spawn.abund[i-n.init.years])
      #Norm Likelihood
      # spawn.abund[i-n.init.years] ~ dnorm(spawn.pred[i-n.init.years], tau.spawn)
      # spawn.abund[i-n.init.years] ~ dlnorm(log(spawn.pred[i-n.init.years]), tau.spawn)
      spawn.abund[i-n.init.years] ~ dnorm(log(spawn.pred[i-n.init.years]), tau.spawn)
      #For lognormal comparison
      mu_spawn.pred[i-n.init.years] <- exp(log(spawn.pred[i-n.init.years]) + sigma.spawn^2/2)
      mode_spawn.pred[i-n.init.years] <- exp(log(spawn.pred[i-n.init.years]) - sigma.spawn^2)
      # sigma_spawn.pred <- sqrt(exp(2*log(spawn.pred)+sigma.spawn^2)*(exp(sigma.spawn^2)-1))
    }#next i
    

    #=======================================================================
    #LIKELIHOOD FOR SPAWNING AGE COMPOSITION
    dirichWeight <- 50#~ dunif(1,100)#50
    # dirichWeight ~ dunif(1,100)#50
#     for(i in 1:2) {
#       temp.dirichWeight[i] ~ dunif(1,100)
#     }
#     sort.dirichWeight[1:2] <- sort(temp.dirichWeight)
#     dirichWeight <- sort.dirichWeight[2]
#     bycatch_dirichWeight <- sort.dirichWeight[1]
    #LIKELIHOOD FOR AGE COMPOSITION
    for(i in (n.init.years+1):(n.years+3)) { #8:31
#       dirichWeight[i] <- total[i]
      for(o in 1:n.oAges) {
        # temp.agecomp.pred[y+stage.offset[o+n.fw.stages],o] <- Esc[y,o]
        # agecomp.pred[i-n.init.years,o] <- Esc.mtx[i,o]/sum(Esc.mtx[i,])
        # adj.Esc.mtx[i,o] <- Esc.mtx[i,o]+0.1 #Works
        # alpha[i-n.init.years,o] <- (adj.Esc.mtx[i,o]/sum(adj.Esc.mtx[i,])) #Works
        agecomp.pred[i-n.init.years,o] <- (Esc.mtx[i,o]/sum(Esc.mtx[i,]))
        # agecomp.obs[i-n.init.years,o] <- agecomp[i,o]
        agecomp.obs[i-n.init.years,o] <- (agecomp[i-n.init.years,o]/sum(agecomp[i-n.init.years,]))
        
        adj.Esc.mtx[i,o] <- Esc.mtx[i,o]#+0.1 #offset for to ensure no dirichlet error is thrown.
        dirich.alpha[i-n.init.years,o] <-  dirichWeight*(adj.Esc.mtx[i,o]/sum(adj.Esc.mtx[i,])+0.001)
        # dirich.alpha[i-n.init.years,o] <-  total[i-n.init.years]*(adj.Esc.mtx[i,o]/sum(adj.Esc.mtx[i,])+0.001)
      }#next o
      #Dirichlet-Multinomial
      # alpha[i-n.init.years,1:n.oAges] ~ ddirich(agecomp.pred[i-n.init.years,1:n.oAges])  #Runs but 
      alpha[i-n.init.years,1:n.oAges] ~ ddirich(dirich.alpha[i-n.init.years,1:n.oAges])#+0.1)
      # alpha[i-n.init.years,1:n.oAges] ~ ddirich(Esc.mtx[i,1:n.oAges]+1) # Fits abund worse
      
      agecomp[i-n.init.years,1:n.oAges] ~ dmulti(alpha[i-n.init.years,1:n.oAges],total[i-n.init.years])
      #Normal Multinomial - Works
      # agecomp[i,1:n.oAges] ~ dmulti(agecomp.pred[i-n.init.years,1:n.oAges],total[i])
    }#next i
    
    #=======================================================================
    #LIKELIHOOD FOR JUVENILE ABUNDANCE DATA
    for(i in 1:n.juvi.locs) {
      juvi.pred[i] <- N[loc.juvi.N[i],4] #4 -those surviving FW stage 2
      #Posterior predictive
      ln.juvi.post.pred[i] ~ dnorm(log(juvi.pred[i]), tau.juvi)
      juvi.post.pred[i] <- exp(ln.juvi.post.pred[i])
      
      # juvi.obs[i] <- juvi.abund[loc.juvi.abund[i]]
      juvi.obs[i] <- exp(juvi.abund[loc.juvi.abund[i]])
      
      # juvi.abund[i] ~ dnorm(juvi.pred[i],tau.juvi)
      # juvi.abund[i] ~ dlnorm(log(juvi.pred[i]),tau.juvi)
      juvi.abund[i] ~ dnorm(log(juvi.pred[i]),tau.juvi)
    }
    
    #=======================================================================
    #LIKELIHOOD FOR BYCATCH ABUNDANCE AND AGE COMP DATA
    #1:20 in Bycatch Data
    #9:28 in bycatch.mtx (y+8)
    bycatch_dirichWeight <- 5#~ dunif(1,100)#5
    # bycatch_dirichWeight ~ dunif(1,100)#5
    # bycatch_dirichWeight <- 5
    for(y in 1:20) { #By calendar year in the bycatch data 
      #Sum of Bycatch by Calendar year
      bycatch.pred[y] <- sum(bycatch.mtx[y+8,])
      #Posterior predictive
      ln.bycatch.post.pred[y] ~ dnorm(log(bycatch.pred[y]), tau.bycatch)
      bycatch.post.pred[y] <- exp(ln.bycatch.post.pred[y])
      # bycatch.obs[y] <- bycatch[y]
      bycatch.obs[y] <- exp(bycatch[y])
      #Extract Observations
      # bycatch.obs[i-n.init.years] <- bycatch[i-n.init.years] #For plotting purpose
      #Norm Likelihood
      # bycatch[y] ~ dlnorm(log(bycatch.pred[y]), tau.bycatch)
      bycatch[y] ~ dnorm(log(bycatch.pred[y]), tau.bycatch)
      
      #Dirichlet-multinomial likelihood for bycatch age composition
      for(o in 1:n.oAges) {
        bycatch_agecomp.pred[y,o] <- (bycatch.mtx[y+8,o]/sum(bycatch.mtx[y+8,]))
        bycatch_agecomp.obs[y,o] <- (ac.bycatch[y,o])/sum(ac.bycatch[y,])
        
        weighted_bycatch_agecomp.pred[y,o] <-  bycatch_dirichWeight*(bycatch_agecomp.pred[y,o]+0.001)
        # weighted_bycatch_agecomp.pred[y,o] <-  sum_ac.bycatch[y]*(bycatch_agecomp.pred[y,o]+0.001)
      }#next o
#       #Dirichlet-Multinomial
      bycatch_alpha[y,1:n.oAges] ~ ddirich(weighted_bycatch_agecomp.pred[y,1:n.oAges])
      ac.bycatch[y,1:n.oAges] ~ dmulti(bycatch_alpha[y,1:n.oAges],sum_ac.bycatch[y])
      #Straight multinomial
      # ac.bycatch[y,1:n.oAges] ~ dmulti(bycatch_agecomp.pred[y,1:n.oAges],sum_ac.bycatch[y])
    }#next i
    

    
  }
  
  parameters.to.save <- c('agecomp.pred', 'agecomp.obs',
                          'spawn.pred','spawn.obs',
                          'juvi.pred', 'juvi.obs',
                          'beta.zero.prod','beta.prod',
                          'beta.zero.cap', #'beta.cap',
                          'exp.beta.zero.cap',
                          'tau.spawn', 'sigma.spawn',
                          'tau.juvi', 'sigma.juvi',
                          # 'N',
                          # 'prod.mtx','cap.mtx',
                          # 'prod', 'cap',
                          # 'Eggs',
                          'surv',
                          # 'Ret',
                          # 'Esc',
                          # 'Catch',
                          'mat.param',
                          'pind',
                          'ind.prod', #'ind.cap'
                          # 'total'
                          # 'Esc.mtx'
                          # 'spawn.pred'
                          # 'temp.mtx')
                          # 'prod.eff','cap.eff'
                          'dirichWeight',
                          
                          'sigma.prod',
                          
                          'sigma.bycatch',
                          'F_bycatch',
                          'bycatch.pred',
                          'bycatch.obs',
                          'mu.Fmort',
                          'sd.Fmort',
                          
                          'bycatch_dirichWeight',
                          'bycatch_agecomp.pred',
                          'bycatch_agecomp.obs',
                          'sel_bycatch',
                          'est.init.spawners',
                          #Derived Quantities of interest for lognormal
                          'mu_spawn.pred',
                          'mode_spawn.pred',
                          # 'sigma_spawn.pred'
                          'slp.beta.zero.prod',
                          'bycatch.post.pred',
                          'spawn.post.pred',
                          'juvi.post.pred'
                          # 'selex_1',
                          # 'selex_2'
)
  
  #Determine which Juvi Brood years to reference
#   juvi.locs <- juvi.abund.years-stage.offset[2]
#   juvi.locs <- juvi.locs[which(juvi.locs<=max(years))]
#   juvi.locs <- which(years %in% juvi.locs)
#   n.juvi.locs <- length(juvi.locs)
  
  
  Data = list('spawn.abund'=log(spawn.abund),
              'juvi.abund'=log(juvi.abund),
              'agecomp'=as.matrix(agecomp),
              'stage.offset'=stage.offset, 'covar.offset'=covar.offset,
              'fecundity'=fecundity, 'femAge'=femAge, 'femProp'=femProp,
              'hr'=hr,
              # 'bycatch'=bycatch,
              'n.stages'=n.stages, 'n.fw.stages'=n.fw.stages, 'n.o.stages'=n.o.stages,
              'n.init.years'=n.init.years, 
              'n.prod.covars'=n.prod.covars, 'n.cap.covars'=n.cap.covars, 'n.covars'=n.covars,
              # 'n.prod.ref'=n.prod.ref, 'n.cap.ref'=n.cap.ref,
              'prod.ref'=prod.ref, 'cap.ref'=cap.ref,
              'covars'=covars,
              'n.years'=n.years, 'n.oAges'=n.oAges, 'n.fem.oAges'=n.fem.oAges,
              'beta.prior.var'=2.5,
              'init.spawners'=rep(6e3,n.init.years),
            
              #Juvenile Data
              'n.juvi.locs'=n.juvi.locs,
              'loc.juvi.N'=loc.juvi.N,
              'loc.juvi.abund'=loc.juvi.abund,
              
              
              'total'=total,
              'dirich.prior'=rep(1,n.oAges),
              'bycatch'=log(bycatch),
              'ac.bycatch'=ac.bycatch,
              'sum_ac.bycatch'=sum_ac.bycatch,
              'meanlog.spawnOE'=meanlog.spawnOE,
              'sdlog.spawnOE'=sdlog.spawnOE
              # 'tau.juvi'=0.1))
  )
  
  # Define initial values
  InitFn = function() {
    # beta.zero.prod <- rep(0,n.stages)
    # beta.zero.prod <- rep(0,n.fw.stages+1)
    # beta.zero.prod <- rnorm(n.fw.stages+1, 0,1)
    # beta.zero.cap <- rep(12,n.fw.stages+1)
    # exp.beta.zero.cap <- runif(2,exp(10),exp(15))
    beta.prod <- rnorm(n.prod.covars,0,1)
    # ind.prod <- rep(1,n.prod.covars)
    ind.prod <- sample(c(0,1), n.prod.covars, replace=TRUE)
    # beta.cap <- rnorm(n.cap.covars,0,1)
    sigma.spawn <- runif(1,1e-2,0.5)
    sigma.juvi <- runif(1,1e-2,0.5)
    sigma.prod <- runif(1,0.1,1)
    sigma.bycatch <- runif(1,1e-2,0.5)
    slp.beta.zero.prod <- 0
    Return <- list(#beta.zero.prod=beta.zero.prod,
                   # exp.beta.zero.cap=exp.beta.zero.cap,
                   beta.prod=beta.prod,
                   # ind.prod=ind.prod
                   # beta.cap=beta.cap
                   sigma.spawn=sigma.spawn,
                   sigma.juvi=sigma.juvi,
                   sigma.prod=sigma.prod,
                   sigma.bycatch=sigma.bycatch,
                   slp.beta.zero.prod=slp.beta.zero.prod,
                   ind.prod=ind.prod
                   )
    return(Return)
  }
  
  #=========================================================================
  ##### RUN JAGS MODEL ######
  
  # out <- NULL
  # out.mcmc <- NULL
  
  #Set Number of Iterations
  # Nsim = Nburnin = sims
  
  Nsim <- sims
  Nburnin <- sims
  out <- FALSE
  
  if(fit.model==TRUE) {
    
    #Single Core
    if(run.parallel==FALSE) {
      
      out <- jags(model.file=yukon_bayes, inits=InitFn, working.directory=NULL, data=Data, parameters.to.save=parameters.to.save,
                  n.chains=chains, n.thin=thins, n.iter=Nsim+Nburnin, n.burnin=Nburnin)   
    }else {
      
      #Parallel
      # system.time(
      # Nsim = Nburnin = sims
      
      #NOT WORK WITHIN FUNCTION - works if run on own
      #   out <- jags.parallel(model.file=yukon_bayes, inits=InitFn, working.directory=NULL, data=Data, parameters.to.save=parameters.to.save,
      #                        n.chains=3, n.thin=thins, n.iter=Nsim+Nburnin, n.burnin=Nburnin,
      #                        export_obj_names=c('thins','Nsim','Nburnin'))   
      
      #NOT WORK WITHIN FUNCTION
      out <- jags.parallel(model.file=yukon_bayes, inits=NULL, working.directory=NULL, data=Data, 
                           parameters.to.save=parameters.to.save,
                           n.chains=chains, n.thin=thins, n.iter=sims+sims, n.burnin=sims,
                           export_obj_names=c('thins','sims','chains'))
      
      #AutoJags
      # out.2 <- autojags(out, n.iter=1000, n.thin=thins)
      # )
    }
    #Save Output
    save(out, file=paste('Output/',pop,'_',model,'.out.RData', sep=''))
    
    #Save speadsheet
    write.csv(out$BUGSoutput$summary, file=paste('Output/',pop,'_',model,'_summary.csv', sep=''))
  }else {
    #Load Output
    load(file=paste('Output/',pop,'_',model,'.out.RData', sep=''))
  }
  
  ##### SECTION FOR UPDATING MCMC CHAIN #####
#   if(update==TRUE) {
      # recompile(out)
    # out.2 <- update(out, n.iter=1e3, n.thin=thins)
#     # save(out, file=paste('Output/',pop,'_',model,'.out.RData', sep=''))
#     # write.csv(out$BUGSoutput$summary, file=paste('Output/',pop,'_',model,'_summary.csv', sep=''))
#   }
  ######## Create mcmc object ########
  out.mcmc <- as.mcmc(out)
  
  #============================================
  pdf(paste('Plots/',pop,'_',model,'_Figs.pdf', sep=''), height=6, width=8)#10
  #PLOT FITS
  par(mfrow=c(2,1), oma=c(0,0,2,0), mar=c(3,4,0,1))
  #Adult
  spawn.plot.yrs <- 1987:2010
  
  y.lim <- c(min(apply(out$BUGSoutput$sims.list$spawn.pred,c(2), quantile, probs=0.025), 
                 apply(out$BUGSoutput$sims.list$spawn.obs,c(2), mean)),
             max(apply(out$BUGSoutput$sims.list$spawn.pred,c(2), quantile, probs=0.975), 
                 apply(out$BUGSoutput$sims.list$spawn.obs,c(2), mean)))
  
  plot(x=spawn.plot.yrs, y=apply(out$BUGSoutput$sims.list$spawn.obs,c(2), mean), type='p', pch=21, bg=rgb(0,0,1, alpha=0.5),
       xlab='', ylab='Spawning Abundance', ylim=y.lim, cex=1.5)
  axis(side=1, at=spawn.plot.yrs, labels=FALSE)
  axis(1, labels=FALSE)
  
  lines(x=spawn.plot.yrs, y=apply(out$BUGSoutput$sims.list$spawn.obs,c(2), mean), col=rgb(0,0,1, alpha=0.5), lwd=2)
  
  #Model
  polygon(x=c(spawn.plot.yrs,rev(spawn.plot.yrs)), y=c(apply(out$BUGSoutput$sims.list$spawn.pred,c(2), quantile, probs=0.975),
                                                       rev(apply(out$BUGSoutput$sims.list$spawn.pred,c(2), quantile, probs=0.025))), 
          col=rgb(1,0,0,alpha=0.25), border=FALSE)
  polygon(x=c(spawn.plot.yrs,rev(spawn.plot.yrs)), y=c(apply(out$BUGSoutput$sims.list$spawn.pred,c(2), quantile, probs=0.75),
                                                       rev(apply(out$BUGSoutput$sims.list$spawn.pred,c(2), quantile, probs=0.25))), 
          col=rgb(1,0,0,alpha=0.25), border=FALSE)
  lines(x=spawn.plot.yrs, y=apply(out$BUGSoutput$sims.list$spawn.pred,c(2), median), col='red', lwd=1)
  
  #Juvi
  juvi.plot.yrs <- juvi.abund.years[1:n.juvi.locs]
  
  y.lim <- c(min(apply(out$BUGSoutput$sims.list$juvi.pred,c(2), quantile, probs=0.025), 
                 apply(out$BUGSoutput$sims.list$juvi.obs,c(2), mean)),
             max(apply(out$BUGSoutput$sims.list$juvi.pred,c(2), quantile, probs=0.975), 
                 apply(out$BUGSoutput$sims.list$juvi.obs,c(2), mean)))
  
  plot(x=juvi.plot.yrs, y=apply(out$BUGSoutput$sims.list$juvi.obs,c(2), mean), type='p', pch=21, bg=rgb(0,0,1, alpha=0.5),
       xlab='', ylab='Survey Index', ylim=y.lim, cex=1.5)
  axis(side=1, at=juvi.plot.yrs, labels=FALSE)
  axis(1, labels=FALSE)
  
  lines(x=juvi.plot.yrs, y=apply(out$BUGSoutput$sims.list$juvi.obs,c(2), mean), col=rgb(0,0,1, alpha=0.5), lwd=2)
  
  #Model
  polygon(x=c(juvi.plot.yrs,rev(juvi.plot.yrs)), y=c(apply(out$BUGSoutput$sims.list$juvi.pred,c(2), quantile, probs=0.975),
                                                     rev(apply(out$BUGSoutput$sims.list$juvi.pred,c(2), quantile, probs=0.025))), 
          col=rgb(1,0,0,alpha=0.25), border=FALSE)
  polygon(x=c(juvi.plot.yrs,rev(juvi.plot.yrs)), y=c(apply(out$BUGSoutput$sims.list$juvi.pred,c(2), quantile, probs=0.75),
                                                     rev(apply(out$BUGSoutput$sims.list$juvi.pred,c(2), quantile, probs=0.25))), 
          col=rgb(1,0,0,alpha=0.25), border=FALSE)
  
  lines(x=juvi.plot.yrs, y=apply(out$BUGSoutput$sims.list$juvi.pred,c(2), median), col='red', lwd=1)
  
  mtext(paste(pop,'River Posterior'), side=3, outer=TRUE, line=0.5)
  
  #========================================================
  #PLOT FIT IN LOG SPACE
  par(mfrow=c(2,1), oma=c(0,0,2,0), mar=c(3,4,0,1))
  #Adult
  spawn.plot.yrs <- 1987:2010
  
  ln.spawn.pred <- log(out$BUGSoutput$sims.list$spawn.pred)
  ln.spawn.obs <- log(out$BUGSoutput$sims.list$spawn.obs)
  
  y.lim <- c(min(apply(ln.spawn.pred,c(2), quantile, probs=0.025), 
                 apply(ln.spawn.obs,c(2), mean)),
             max(apply(ln.spawn.pred,c(2), quantile, probs=0.975), 
                 apply(ln.spawn.obs,c(2), mean)))
  
  plot(x=spawn.plot.yrs, y=apply(ln.spawn.obs,c(2), mean), type='p', pch=21, bg=rgb(0,0,1, alpha=0.5),
       xlab='', ylab='ln(Spawning Abundance)', ylim=y.lim, cex=1.5)
  axis(side=1, at=spawn.plot.yrs, labels=FALSE)
  axis(1, labels=FALSE)
  
  lines(x=spawn.plot.yrs, y=apply(ln.spawn.obs,c(2), mean), col=rgb(0,0,1, alpha=0.5), lwd=2)
  
  #Model
  polygon(x=c(spawn.plot.yrs,rev(spawn.plot.yrs)), y=c(apply(ln.spawn.pred,c(2), quantile, probs=0.975),
                                                       rev(apply(ln.spawn.pred,c(2), quantile, probs=0.025))), 
          col=rgb(1,0,0,alpha=0.25), border=FALSE)
  
  polygon(x=c(spawn.plot.yrs,rev(spawn.plot.yrs)), y=c(apply(ln.spawn.pred,c(2), quantile, probs=0.75),
                                                       rev(apply(ln.spawn.pred,c(2), quantile, probs=0.25))), 
          col=rgb(1,0,0,alpha=0.25), border=FALSE)
  lines(x=spawn.plot.yrs, y=apply(ln.spawn.pred,c(2), median), col='red', lwd=1)
  
  #Juvi
  juvi.plot.yrs <- juvi.abund.years[1:n.juvi.locs]
  
  ln.juvi.pred <- log(out$BUGSoutput$sims.list$juvi.pred)
  ln.juvi.obs <- log(out$BUGSoutput$sims.list$juvi.obs)
  
  y.lim <- c(min(apply(ln.juvi.pred,c(2), quantile, probs=0.025), 
                 apply(ln.juvi.obs,c(2), mean)),
             max(apply(ln.juvi.pred,c(2), quantile, probs=0.975), 
                 apply(ln.juvi.obs,c(2), mean)))
  
  plot(x=juvi.plot.yrs, y=apply(ln.juvi.obs,c(2), mean), type='p', pch=21, bg=rgb(0,0,1, alpha=0.5),
       xlab='', ylab='ln(Survey Index)', ylim=y.lim, cex=1.5)
  axis(side=1, at=juvi.plot.yrs, labels=FALSE)
  axis(1, labels=FALSE)
  
  lines(x=juvi.plot.yrs, y=apply(ln.juvi.obs,c(2), mean), col=rgb(0,0,1, alpha=0.5), lwd=2)
  
  #Model
  polygon(x=c(juvi.plot.yrs,rev(juvi.plot.yrs)), y=c(apply(ln.juvi.pred,c(2), quantile, probs=0.975),
                                                     rev(apply(ln.juvi.pred,c(2), quantile, probs=0.025))), 
          col=rgb(1,0,0,alpha=0.25), border=FALSE)
  polygon(x=c(juvi.plot.yrs,rev(juvi.plot.yrs)), y=c(apply(ln.juvi.pred,c(2), quantile, probs=0.75),
                                                     rev(apply(ln.juvi.pred,c(2), quantile, probs=0.25))), 
          col=rgb(1,0,0,alpha=0.25), border=FALSE)
  
  lines(x=juvi.plot.yrs, y=apply(ln.juvi.pred,c(2), median), col='red', lwd=1)
  mtext(paste(pop,'River Posterior'), side=3, outer=TRUE, line=0.5)
  
  ##############################################################################################################
  #POSTERIOR PREDICTIVE DISTRIBUTION
  par(mfrow=c(2,1), oma=c(0,0,2,0), mar=c(3,4,0,1))
  #Adult
  spawn.plot.yrs <- 1987:2010
  
  y.lim <- c(min(apply(out$BUGSoutput$sims.list$spawn.post.pred,c(2), quantile, probs=0.025), 
                 apply(out$BUGSoutput$sims.list$spawn.obs,c(2), mean)),
             max(apply(out$BUGSoutput$sims.list$spawn.post.pred,c(2), quantile, probs=0.975), 
                 apply(out$BUGSoutput$sims.list$spawn.obs,c(2), mean)))
  
  plot(x=spawn.plot.yrs, y=apply(out$BUGSoutput$sims.list$spawn.obs,c(2), mean), type='p', pch=21, bg=rgb(0,0,1, alpha=0.5),
       xlab='', ylab='Spawning Abundance', ylim=y.lim, cex=1.5)
  axis(side=1, at=spawn.plot.yrs, labels=FALSE)
  axis(1, labels=FALSE)
  
  lines(x=spawn.plot.yrs, y=apply(out$BUGSoutput$sims.list$spawn.obs,c(2), mean), col=rgb(0,0,1, alpha=0.5), lwd=2)
  
  #Model
  polygon(x=c(spawn.plot.yrs,rev(spawn.plot.yrs)), y=c(apply(out$BUGSoutput$sims.list$spawn.post.pred,c(2), quantile, probs=0.975),
                                                       rev(apply(out$BUGSoutput$sims.list$spawn.post.pred,c(2), quantile, probs=0.025))), 
          col=rgb(1,0,0,alpha=0.25), border=FALSE)
  polygon(x=c(spawn.plot.yrs,rev(spawn.plot.yrs)), y=c(apply(out$BUGSoutput$sims.list$spawn.post.pred,c(2), quantile, probs=0.75),
                                                       rev(apply(out$BUGSoutput$sims.list$spawn.post.pred,c(2), quantile, probs=0.25))), 
          col=rgb(1,0,0,alpha=0.25), border=FALSE)
  lines(x=spawn.plot.yrs, y=apply(out$BUGSoutput$sims.list$spawn.post.pred,c(2), median), col='red', lwd=1)
  
  #Juvi
  juvi.plot.yrs <- juvi.abund.years[1:n.juvi.locs]
  
  y.lim <- c(min(apply(out$BUGSoutput$sims.list$juvi.post.pred,c(2), quantile, probs=0.025), 
                 apply(out$BUGSoutput$sims.list$juvi.obs,c(2), mean)),
             max(apply(out$BUGSoutput$sims.list$juvi.post.pred,c(2), quantile, probs=0.975), 
                 apply(out$BUGSoutput$sims.list$juvi.obs,c(2), mean)))
  
  plot(x=juvi.plot.yrs, y=apply(out$BUGSoutput$sims.list$juvi.obs,c(2), mean), type='p', pch=21, bg=rgb(0,0,1, alpha=0.5),
       xlab='', ylab='Survey Index', ylim=y.lim, cex=1.5)
  axis(side=1, at=juvi.plot.yrs, labels=FALSE)
  axis(1, labels=FALSE)
  
  lines(x=juvi.plot.yrs, y=apply(out$BUGSoutput$sims.list$juvi.obs,c(2), mean), col=rgb(0,0,1, alpha=0.5), lwd=2)
  
  #Model
  polygon(x=c(juvi.plot.yrs,rev(juvi.plot.yrs)), y=c(apply(out$BUGSoutput$sims.list$juvi.post.pred,c(2), quantile, probs=0.975),
                                                     rev(apply(out$BUGSoutput$sims.list$juvi.post.pred,c(2), quantile, probs=0.025))), 
          col=rgb(1,0,0,alpha=0.25), border=FALSE)
  polygon(x=c(juvi.plot.yrs,rev(juvi.plot.yrs)), y=c(apply(out$BUGSoutput$sims.list$juvi.post.pred,c(2), quantile, probs=0.75),
                                                     rev(apply(out$BUGSoutput$sims.list$juvi.post.pred,c(2), quantile, probs=0.25))), 
          col=rgb(1,0,0,alpha=0.25), border=FALSE)
  
  lines(x=juvi.plot.yrs, y=apply(out$BUGSoutput$sims.list$juvi.post.pred,c(2), median), col='red', lwd=1)
  
  mtext(paste(pop,'River Posterior Predictive Distribution'), side=3, outer=TRUE, line=0.5)
  
  #========================================================
  #PLOT FIT IN LOG SPACE
  par(mfrow=c(2,1), oma=c(0,0,2,0), mar=c(3,4,0,1))
  #Adult
  spawn.plot.yrs <- 1987:2010
  
  ln.spawn.post.pred <- log(out$BUGSoutput$sims.list$spawn.post.pred)
  ln.spawn.obs <- log(out$BUGSoutput$sims.list$spawn.obs)
  
  y.lim <- c(min(apply(ln.spawn.post.pred,c(2), quantile, probs=0.025), 
                 apply(ln.spawn.obs,c(2), mean)),
             max(apply(ln.spawn.post.pred,c(2), quantile, probs=0.975), 
                 apply(ln.spawn.obs,c(2), mean)))
  
  plot(x=spawn.plot.yrs, y=apply(ln.spawn.obs,c(2), mean), type='p', pch=21, bg=rgb(0,0,1, alpha=0.5),
       xlab='', ylab='ln(Spawning Abundance)', ylim=y.lim, cex=1.5)
  axis(side=1, at=spawn.plot.yrs, labels=FALSE)
  axis(1, labels=FALSE)
  
  lines(x=spawn.plot.yrs, y=apply(ln.spawn.obs,c(2), mean), col=rgb(0,0,1, alpha=0.5), lwd=2)
  
  #Model
  polygon(x=c(spawn.plot.yrs,rev(spawn.plot.yrs)), y=c(apply(ln.spawn.post.pred,c(2), quantile, probs=0.975),
                                                       rev(apply(ln.spawn.post.pred,c(2), quantile, probs=0.025))), 
          col=rgb(1,0,0,alpha=0.25), border=FALSE)
  
  polygon(x=c(spawn.plot.yrs,rev(spawn.plot.yrs)), y=c(apply(ln.spawn.post.pred,c(2), quantile, probs=0.75),
                                                       rev(apply(ln.spawn.post.pred,c(2), quantile, probs=0.25))), 
          col=rgb(1,0,0,alpha=0.25), border=FALSE)
  lines(x=spawn.plot.yrs, y=apply(ln.spawn.post.pred,c(2), median), col='red', lwd=1)
  
  #Juvi
  juvi.plot.yrs <- juvi.abund.years[1:n.juvi.locs]
  
  ln.juvi.post.pred <- log(out$BUGSoutput$sims.list$juvi.post.pred)
  ln.juvi.obs <- log(out$BUGSoutput$sims.list$juvi.obs)
  
  y.lim <- c(min(apply(ln.juvi.post.pred,c(2), quantile, probs=0.025), 
                 apply(ln.juvi.obs,c(2), mean)),
             max(apply(ln.juvi.post.pred,c(2), quantile, probs=0.975), 
                 apply(ln.juvi.obs,c(2), mean)))
  
  plot(x=juvi.plot.yrs, y=apply(ln.juvi.obs,c(2), mean), type='p', pch=21, bg=rgb(0,0,1, alpha=0.5),
       xlab='', ylab='ln(Survey Index)', ylim=y.lim, cex=1.5)
  axis(side=1, at=juvi.plot.yrs, labels=FALSE)
  axis(1, labels=FALSE)
  
  lines(x=juvi.plot.yrs, y=apply(ln.juvi.obs,c(2), mean), col=rgb(0,0,1, alpha=0.5), lwd=2)
  
  #Model
  polygon(x=c(juvi.plot.yrs,rev(juvi.plot.yrs)), y=c(apply(ln.juvi.post.pred,c(2), quantile, probs=0.975),
                                                     rev(apply(ln.juvi.post.pred,c(2), quantile, probs=0.025))), 
          col=rgb(1,0,0,alpha=0.25), border=FALSE)
  polygon(x=c(juvi.plot.yrs,rev(juvi.plot.yrs)), y=c(apply(ln.juvi.post.pred,c(2), quantile, probs=0.75),
                                                     rev(apply(ln.juvi.post.pred,c(2), quantile, probs=0.25))), 
          col=rgb(1,0,0,alpha=0.25), border=FALSE)
  
  lines(x=juvi.plot.yrs, y=apply(ln.juvi.post.pred,c(2), median), col='red', lwd=1)
  
  mtext(paste(pop,'River Posterior Predictive Distribution'), side=3, outer=TRUE, line=0.5)
  ##############################################################################################################
  
  #=================================================================================
  #AGE COMPOSITION FITS
  cols <- rainbow(n=n.oAges)
  par(mfrow=c(2,1), oma=c(2,2,2,1), mar=c(3,4,1,1))
  ac.obs <- apply(out$BUGSoutput$sims.list$agecomp.obs,c(2,3), mean)
  prop.obs <- ac.obs/apply(ac.obs, 1, sum)
  barplot(t(prop.obs), col=cols, ylab='Observed', las=2)
  
  ac.pred <- apply(out$BUGSoutput$sims.list$agecomp.pred,c(2,3), median)
  prop.pred <- ac.pred/apply(ac.pred,1,sum)
  
  barplot(t(prop.pred), col=cols, ylab='Predicted', las=2, names.arg=1987:(1987+nrow(prop.pred)-1))
  
  mtext('Age Composition Proportion', side=2, outer=TRUE, font=2, cex=1.25, line=0.5)
  mtext('Year', side=1, outer=TRUE, font=2, cex=1.25, line=0.5)
  
  par(mfrow=c(1,1), mar=c(4,4,1,1), oma=c(0,0,0,0))
  
  lim <- c(min(prop.pred, prop.obs), max(prop.pred, prop.obs))
  
  plot(x=NULL, y=NULL, xlab='Observed Proportion', ylab='Predicted Proportion', main=pop, xlim=lim, ylim=lim)
  a <- 1
  for(a in 1:n.oAges) {
    points(x=prop.obs[,a], y=prop.pred[,a], pch=21, bg=cols[a])
  }#next a
  segments(x0=-1, y0=-1, x1=2, y1=2, lty=2, lwd=2)
  legend('topleft', legend=oAges, pch=21, pt.bg=cols, title='Ocean Age')
  #=============================================
  # out.surv <- apply(out$BUGSoutput$sims.list$surv,c(2,3), median)
  # out.prod <- apply(out$BUGSoutput$sims.list$prod,c(2,3), median)
  # out.cap <- apply(out$BUGSoutput$sims.list$cap,c(2,3), median)
  # 
  # out.Esc <- apply(out$BUGSoutput$sims.list$Esc,c(2,3), median)
  # out.Ret <- apply(out$BUGSoutput$sims.list$Ret,c(2,3), median)
  # out.Catch <- apply(out$BUGSoutput$sims.list$Catch,c(2,3), median)
  # 
  # out.agecomp.pred <- apply(out$BUGSoutput$sims.list$agecomp.pred,c(2,3), median)
  # out.agecomp.obs <- apply(out$BUGSoutput$sims.list$agecomp.obs,c(2,3), median)
  # 
  # out.temp.agecomp.pred <- apply(out$BUGSoutput$sims.list$temp.agecomp.pred,c(2,3), median)
  
  #Plot Survival
  par(mfrow=c(n.stages,1), oma=c(3,3,2,1), mar=c(0,4,0,0))
  temp.stage.names <- stage.names
  temp.stage.names[2] <- 'Downstream\nNearshore'
  s <- 1
  for(s in 1:n.stages) {
    plot(x=years, y=apply(out$BUGSoutput$sims.list$surv,c(2,3), median)[,s], ylab=temp.stage.names[s], 
         type='l', xaxt='n', ylim=c(0,1), col='red', lwd=2)
    polygon(x=c(years,rev(years)), y=c(apply(out$BUGSoutput$sims.list$surv,c(2,3), quantile, probs=0.025)[,s],
                                       rev(apply(out$BUGSoutput$sims.list$surv,c(2,3), quantile, probs=0.975)[,s])),
            col=rgb(1,0,0, alpha=0.25), border=FALSE)
    
    polygon(x=c(years,rev(years)), y=c(apply(out$BUGSoutput$sims.list$surv,c(2,3), quantile, probs=0.25)[,s],
                                       rev(apply(out$BUGSoutput$sims.list$surv,c(2,3), quantile, probs=0.75)[,s])),
            col=rgb(1,0,0, alpha=0.25), border=FALSE)
  }#next s
  axis(1, at=years, labels=years)
  mtext('Realized Survival', side=2, outer=TRUE, line=2)
  
  #=============================================
  par(mfrow=c(1,1), mar=c(4,7,1,2), oma=c(0,0,0,0))
  #Beta Prod
  prod.labs <- paste(prod.covar.names, '\np=', round(apply(out$BUGSoutput$sims.list$ind.prod, 2, mean),2), sep='')
  caterplot(out.mcmc, parms='beta.prod', collapse=TRUE, reorder=FALSE, 
            quantiles=list(outer=c(0.025,0.975),inner=c(0.25,0.75)), labels=prod.labs)
  mtext('Beta Prod', side=1, line=2.5)
  abline(v=0, col='red')
  
  prod.labs <- paste(prod.covar.names, '\np=', round(apply(out$BUGSoutput$sims.list$ind.prod, 2, mean),2), sep='')
  caterplot(out.mcmc, parms='beta.prod', collapse=FALSE, reorder=FALSE, 
            quantiles=list(outer=c(0.025,0.975),inner=c(0.25,0.75)), labels=prod.labs)
  mtext('Beta Prod', side=1, line=2.5)
  abline(v=0, col='red')
  
  #Traceplot indprod
  
  
  #=============================================
  #GGPLOT - Beta prod
  list.beta.prod <- melt(out$BUGSoutput$sims.list$beta.prod) 
  names(list.beta.prod) <- c('sim','var.num','value')
  var.name <- prod.covar.names[list.beta.prod$var.num]
  var.stage <- covar.stage[list.beta.prod$var.num]
  
  list.beta.prod <- data.frame(list.beta.prod, var.name, var.stage)
  list.beta.prod$var.name <- factor(list.beta.prod$var.name, levels=rev(prod.covar.names), ordered=TRUE)
  list.beta.prod$var.stage <- factor(list.beta.prod$var.stage,levels=stage.names, ordered=TRUE)
  
  #Determine which samples to remove p(b|ind.prod==0)
  list.ind.prod <- melt(out$BUGSoutput$sims.list$ind.prod)
  loc.remove <- which(list.ind.prod$value==0)
  list.beta.prod <- list.beta.prod[-loc.remove,]
  
  #Plot it Out
  tmp <- ggplot(list.beta.prod, aes(var.name, value, fill=var.stage))
  tmp <- tmp + ylab('Effect') + xlab('Covariate')
  if(pop=='chena') {
    tmp <- tmp + ggtitle('Chena River')	
  }else {
    tmp <- tmp + ggtitle('Salcha')	
  }
  tmp <- tmp + geom_hline(yintercept=0, col='red')
  tmp <- tmp + scale_y_continuous(limits = quantile(list.beta.prod$value, c(0.001, 0.999)))
  # tmp <- tmp + geom_boxplot(width=0.25, lwd=0.5, outlier.shape=NA)#width=0.1, lwd=0.1
  tmp <- tmp + geom_boxplot(lwd=0.5, outlier.shape=NA)
  tmp <- tmp + coord_flip()
  plot(tmp)
  tmp <- tmp + theme(legend.position='none')
  plot(tmp)
  tmp <- tmp + geom_boxplot(width=0.25, lwd=0.5, outlier.shape=NA)
  tmp <- tmp + geom_violin(alpha = 0.5, lwd=0.1, scale='width') 
  plot(tmp)
  
  #With Violins
  tmp <- ggplot(list.beta.prod, aes(var.name, value, fill=var.stage))
  tmp <- tmp + ylab('Effect') + xlab('Covariate')
  tmp <- tmp + ylab('Effect') + xlab('Covariate')
  if(pop=='chena') {
    tmp <- tmp + ggtitle('Chena River')	
  }else {
    tmp <- tmp + ggtitle('Salcha')	
  }
  tmp <- tmp + geom_hline(yintercept=0, col='red')
  tmp <- tmp + scale_y_continuous(limits = quantile(list.beta.prod$value, c(0.001, 0.999)))
  # tmp <- tmp + geom_boxplot(width=0.25, lwd=0.5, outlier.shape=NA)#width=0.1, lwd=0.1
  # tmp <- tmp + geom_boxplot(lwd=0.5, outlier.shape=NA)
  tmp <- tmp + geom_boxplot(width=0.25, lwd=0.5, outlier.shape=NA)
  tmp <- tmp + geom_violin(alpha = 0.5, lwd=0.1, scale='width') 
  tmp <- tmp + coord_flip()
  plot(tmp)
  tmp <- tmp + theme(legend.position='none')
  plot(tmp)

  
  #Separate by stage
  tmp <- ggplot(list.beta.prod, aes(var.name, value, fill=var.stage))
  tmp <- tmp + ylab('Effect') + xlab('Covariate')
  tmp <- tmp + theme(legend.position='none')
  tmp <- tmp + facet_wrap(~var.stage, ncol=1, shrink=TRUE, drop=TRUE)
  # tmp <- tmp + facet_grid(var.stage ~ ., drop=TRUE)
  # tmp <- tmp + facet_grid(var.stage~., shrink=TRUE, drop=TRUE)
  if(pop=='chena') {
    tmp <- tmp + ggtitle('Chena River')	
  }else {
    tmp <- tmp + ggtitle('Salcha')	
  }
  tmp <- tmp + geom_hline(yintercept=0, col='red')
  # tmp <- tmp + geom_boxplot(width=0.25, lwd=0.5, outlier.shape=NA)#width=0.1, lwd=0.1
  tmp <- tmp + geom_boxplot(lwd=0.5, outlier.shape=NA)
  tmp <- tmp + coord_flip()
  plot(tmp)
  tmp <- tmp + geom_violin(alpha = 0.5, lwd=0.1, scale='width') 
  plot(tmp)
  #Separate by s
  
  
  #Beta Cap
#   cap.labs <- paste(cap.covar.names, '\np=', round(apply(out$BUGSoutput$sims.list$ind.cap, 2, mean),2), sep='')
#   caterplot(out.mcmc, parms='beta.cap', collapse=TRUE, reorder=FALSE, 
#             quantiles=list(outer=c(0.025,0.975),inner=c(0.25,0.75)), labels=cap.labs)
#   mtext('Beta Cap', side=1, line=2.5)
#   abline(v=0, col='red')
#   
#   cap.labs <- paste(cap.covar.names, '\np=', round(apply(out$BUGSoutput$sims.list$ind.cap, 2, mean),2), sep='')
#   caterplot(out.mcmc, parms='beta.cap', collapse=FALSE, reorder=FALSE, 
#             quantiles=list(outer=c(0.025,0.975),inner=c(0.25,0.75)), labels=cap.labs)
#   mtext('Beta Cap', side=1, line=2.5)
#   abline(v=0, col='red')
#   
  #=============================================
  # out.Esc.mtx <- apply(out$BUGSoutput$sims.list$Esc.mtx,c(2,3), median)
  
  #Auto Update
  # out.auto <- autojags(out.auto)
  par(mfrow=c(1,1), mar=c(4,6,1,2), oma=c(0,0,0,0))
  
  #Beta Zero Prod
  caterplot(out.mcmc, parms='beta.zero.prod', collapse=TRUE, reorder=FALSE, 
            quantiles=list(outer=c(0.025,0.975),inner=c(0.25,0.75)), labels=temp.stage.names)
  mtext('Beta Zero Prod', side=1, line=2.5)
  
  caterplot(out.mcmc, parms='beta.zero.prod', collapse=FALSE, reorder=FALSE, 
            quantiles=list(outer=c(0.025,0.975),inner=c(0.25,0.75)), labels=temp.stage.names)
  mtext('Beta Zero Prod', side=1, line=2.5)
  
  #Beta Zero Cap
  caterplot(out.mcmc, parms='beta.zero.cap', collapse=TRUE, reorder=FALSE, 
            quantiles=list(outer=c(0.025,0.975),inner=c(0.25,0.75)))#, labels=temp.stage.names[1:(n.fw.stages+1)])
  mtext('Beta Zero Cap', side=1, line=2.5)
  
  caterplot(out.mcmc, parms='exp.beta.zero.cap', collapse=FALSE, reorder=FALSE, 
            quantiles=list(outer=c(0.025,0.975),inner=c(0.25,0.75)))#, labels=temp.stage.names[1:(n.fw.stages+1)])
  mtext('Exp Beta Zero Cap', side=1, line=2.5)
  
  #Reordered
  caterplot(out.mcmc, parms='beta.prod', collapse=TRUE, reorder=TRUE, 
            quantiles=list(outer=c(0.025,0.975),inner=c(0.25,0.75)))#, labels=prod.covar.names)
  abline(v=0, col='red')
#   caterplot(out.mcmc, parms='beta.cap', collapse=TRUE, reorder=TRUE, 
#             quantiles=list(outer=c(0.025,0.975),inner=c(0.25,0.75)), labels=cap.covar.names)
#   abline(v=0, col='red')
  
  
  #===========================================================
  #PLOT INCLUSION PROBABILITY
  par(mfrow=c(1,1), mar=c(4,9,2,1), oma=c(0,0,0,0))
  probs.prod <- apply(out$BUGSoutput$sims.list$ind.prod, c(2), mean)
  # probs.cap <- apply(out$BUGSoutput$sims.list$ind.cap, c(2), mean)
  
  #Prob
  barplot(probs.prod, horiz=TRUE, names.arg=prod.covar.names, las=2, main=pop, col=rgb(0,0,1, alpha=0.5), 
          border='black', xlab='Productivity Inclusion Probability', xlim=c(0,1))
  # abline(v=0, lwd=2)
  abline(v=pretty(c(0,1)), lty=2)
  box(lwd=2)
  abline(v=median(out$BUGSoutput$sims.list$pind), col=rgb(1,0,0,alpha=0.5), lwd=3)
  #Using ggplot
  
  list.ind.prod <- melt(out$BUGSoutput$sims.list$ind.prod)
  names(list.ind.prod) <- c('sim','var.num','value')
  var.name <- prod.covar.names[list.ind.prod$var.num]
  var.stage <- covar.stage[list.ind.prod$var.num]
  
  list.ind.prod <- data.frame(list.ind.prod, var.name, var.stage)
  list.ind.prod$var.name <- factor(list.ind.prod$var.name, levels=rev(prod.covar.names), ordered=TRUE)
  list.ind.prod$var.stage <- factor(list.ind.prod$var.stage,levels=stage.names, ordered=TRUE)
  #Add the pind median to list
  med.pind <- median(out$BUGSoutput$sims.list$pind)
  list.ind.prod <- data.frame(list.ind.prod,med.pind)
  
  #Plot it Out
  tmp <- ggplot(list.ind.prod, aes(var.name, fill=var.stage, alpha=0.75))
  tmp <- tmp + stat_summary_bin(aes(y = value), fun.y = "mean", geom = "bar")
  tmp <- tmp + ylim(0,1)
  tmp <- tmp + coord_flip()
  tmp <- tmp + ylab('Inclusion Probability') + xlab('Covariate')

  # tmp <- tmp + scale_y_continuous(c(0,1))
  if(pop=='chena') {
    tmp <- tmp + ggtitle('Chena River')	
  }else {
    tmp <- tmp + ggtitle('Salcha')	
  }
  
  #Add In prior expectation
  
  plot(tmp)
  tmp <- tmp + theme(legend.position='none')
  plot(tmp)
  tmp <- tmp + geom_abline(intercept = med.pind, slope=0, alpha=0.5, color='red')
  plot(tmp)
  
  #Traceplot for inclusion probability
  traplot(out.mcmc, parms='ind.prod')
  
  #Cap
#   barplot(probs.cap, horiz=TRUE, names.arg=cap.covar.names, las=2, main=pop, col=rgb(0,0,1, alpha=0.5), 
#           border='black', xlab='Capacity Inclusion Probability', xlim=c(0,1))
#   abline(v=0, lwd=2)
  #   apply(out$BUGSoutput$sims.list$ind.prod,c(2), mean)
  #   apply(out$BUGSoutput$sims.list$ind.cap,c(2), mean)
  
  denplot(out.mcmc, parms='beta.prod')
  denplot(out.mcmc, parms='mat.param')
  
  denplot(out.mcmc, parms='beta.zero.prod')
  denplot(out.mcmc, parms='slp.beta.zero.prod')
  
  par(mfrow=c(3,3), oma=c(0,0,2,0), mar=c(5,2,1,2))
  for(i in 1:n.stages) {
    plotPost(inv.logit(out$BUGSoutput$sims.list$beta.zero.prod[,i]), xlab=stage.names[i], xlim=c(0,1))
  }
  mtext('BASELINE SURVIVAL: inv.logit(Beta Zero Prod)', outer=TRUE, line=0)
  
  denplot(out.mcmc, parms='beta.zero.cap') #Edited
  denplot(out.mcmc, parms='exp.beta.zero.cap')
  
  denplot(out.mcmc, parms=list('pind'))
  
  denplot(out.mcmc, parms='sigma.spawn')
  denplot(out.mcmc, parms='sigma.juvi')
  
  # caterplot(out.mcmc, parms='ind.prod')
  
  #========================================================
  #Plot Maturation Schedule
  par(mfrow=c(1,1), oma=c(0,0,0,0), mar=c(4,4,2,1))
  temp.matPar <- as.vector(out$BUGSoutput$sims.list$mat.param)
  temp.ages <- seq(from=1, to=5, length.out=100)
  temp.mat <- matrix(nrow=length(temp.matPar), ncol=length(temp.ages))
  
  i <- 1
  for(i in 1:length(temp.ages)) {
    temp.mat[,i] <-  1/(1+exp(temp.matPar*(5-temp.ages[i])+log(1/0.99-1)))
  }
  plot(x=NULL, y=NULL, type='l', ylim=c(0,1), xlim=c(min(temp.ages), max(temp.ages)),
         main=pop, xlab='Ocean Age', ylab='Maturation Probability')
  polygon(x=c(temp.ages,rev(temp.ages)), y=c(apply(temp.mat, c(2), quantile, probs=0.975),
                                             rev(apply(temp.mat, c(2), quantile, probs=0.025))),
          col=rgb(0,0,1, alpha=0.25), border=FALSE)
  polygon(x=c(temp.ages,rev(temp.ages)), y=c(apply(temp.mat, c(2), quantile, probs=0.75),
                                             rev(apply(temp.mat, c(2), quantile, probs=0.25))),
          col=rgb(0,0,1, alpha=0.25), border=FALSE)
  lines(x=temp.ages, y=apply(temp.mat, c(2), quantile, probs=0.5), lwd=2)
  grid(col='black')
  
  #========================================================
  #Realized Survival
  surv.list <- melt(out$BUGSoutput$sims.list$surv)
  names(surv.list) <- c('n','year','stage.num','value')
  stage <- vector(length=nrow(surv.list))
  temp.stage.names <- stage.names
  temp.stage.names[2] <- 'Downstream\nNearshore'
  stage <- temp.stage.names[surv.list$stage]
  surv.list.final <- data.frame(surv.list,stage)
  surv.list.final$stage <- factor(surv.list.final$stage, levels=temp.stage.names, ordered=TRUE)

  tmp <- ggplot(surv.list.final, aes(stage, value, fill=stage))
  # tmp <- tmp + geom_violin(alpha = 0.5, lwd=0.1, scale='width') 
  tmp <- tmp + ylab('Realized Survival') + xlab('Stage')
  if(pop=='chena') {
    tmp <- tmp + ggtitle('Chena River')	
  }else {
    tmp <- tmp + ggtitle('Salcha')	
  }
  tmp <- tmp + geom_boxplot(width=0.25, lwd=0.5, outlier.shape=NA)#width=0.1, lwd=0.1
  plot(tmp)
  tmp <- tmp + theme(legend.position='none')
  plot(tmp)

  #Plot Average Survival
  avg.surv.list <- melt(apply(out$BUGSoutput$sims.list$surv, c(1,3), mean))
  names(avg.surv.list) <- c('n','stage.num','value')
  stage <- vector(length=nrow(avg.surv.list))
  temp.stage.names <- stage.names
  temp.stage.names[2] <- 'Downstream\nNearshore'
  stage <- temp.stage.names[avg.surv.list$stage]
  avg.surv.list <- data.frame(avg.surv.list,stage)
  avg.surv.list$stage <- factor(avg.surv.list$stage, levels=temp.stage.names, ordered=TRUE)
  
  tmp <- ggplot(avg.surv.list, aes(stage, value, fill=stage))
  # tmp <- tmp + geom_violin(alpha = 0.5, lwd=0.1, scale='width') 
  tmp <- tmp + ylab('Average Survival') + xlab('Stage')
  if(pop=='chena') {
    tmp <- tmp + ggtitle('Chena River')	
  }else {
    tmp <- tmp + ggtitle('Salcha')	
  }
  tmp <- tmp + geom_boxplot(width=0.25, lwd=0.5, outlier.shape=NA)#width=0.1, lwd=0.1
  plot(tmp)
  tmp <- tmp + theme(legend.position='none')
  plot(tmp)
  
  #===========================================================
  #PLOT REALIZED SURVIVAL DIFFERENCES - SIMPLE
  #Total number of posterior samples
  n.samp <- nrow(out$BUGSoutput$sims.list$ind.prod)
  
  #Gatther posterior samples
  out$parameters.to.save
  
  sim.beta.zero.prod <- out$BUGSoutput$sims.list$beta.zero.prod
  sim.beta.zero.cap <- out$BUGSoutput$sims.list$beta.zero.cap
  sim.beta.prod <- out$BUGSoutput$sims.list$beta.prod
  sim.beta.cap <- NULL
  sim.ind.prod <- out$BUGSoutput$sims.list$ind.prod
  sim.mat.param  <- out$BUGSoutput$sims.list$mat.param
  sim.hr <- mean(hr)
  #Simulated Bycatch Mortality
  sim.fmort <- out$BUGSoutput$sims.list$mu.Fmort
  sim.sel_bycatch <- out$BUGSoutput$sims.list$sel_bycatch
  
  #Calculate initial conditions
  sim.spawners <- mean(spawn.abund, na.rm=TRUE)
  avg.femProp <- mean(femProp)
  avg.fecundity <- mean(fecundity)
  avg.femAge <- apply(femAge, 2, mean)
  avg.femAge <- avg.femAge/sum(avg.femAge) #Standardize
  #Calculate initial conditions
  sim.eggs <- sim.spawners*avg.femProp*sum(avg.fecundity*avg.femAge)
  
  #Set up output objects
  sim.N <- array(dim=c(n.stages+2, n.samp, n.covars), dimnames=list(c('Spawners','Eggs',stage.names), c(1:n.samp), covar.names))
  base.N <- array(dim=c(n.stages+2, n.samp), dimnames=list(c('Spawners','Eggs',stage.names), c(1:n.samp)))
  
  #No Comm Fish
  sim.rec <- array(data=0, dim=c(n.samp, n.covars), dimnames=list(c(1:n.samp), covar.names))
  base.rec <- rep(0,n.samp)
  
  #With Comm Fish
  sim.ret <- array(data=0, dim=c(n.samp, n.covars), dimnames=list(c(1:n.samp), covar.names))
  base.ret <- rep(0,n.samp)

  i <- 1
  for(i in 1:n.samp) {  
    c <- 1
    for(c in 1:n.covars) {
      #Create temporary covariate set
      temp.covars <- rep(0,n.covars)
      temp.covars[c] <- 1
      #===============================
      #COVARS
      
      #Initialize Spawners
      sim.N[1,i,c] <- sim.spawners
      #Initialize Eggs
      sim.N[2,i,c] <- sim.eggs

      #Freshwater
      s <- 1
      for(s in 1:n.fw.stages) {
        #Productivity
        temp.prod.loc <- prod.ref[s,which(prod.ref[s,]!=0)]
        if(length(temp.prod.loc)>0) {
          temp.prod <- 1/(1+exp(-1*sim.beta.zero.prod[i,s] - sum(sim.ind.prod[i,temp.prod.loc] * sim.beta.prod[i,temp.prod.loc] * temp.covars[temp.prod.loc])))
          # temp.prod <- 1/(1+exp(-1*sim.beta.zero.prod[i,s] - sum(sim.beta.prod[i,temp.prod.loc] * temp.covars[temp.prod.loc])))
        }else {
          temp.prod <- 1/(1+exp(-1*sim.beta.zero.prod[i,s]))
        }
        #Capacity
        temp.cap.loc <- cap.ref[s,which(cap.ref[s,]!=0)]
        if(length(temp.cap.loc)>0) {
          temp.cap <- exp(sim.beta.zero.cap[i,1] + sum(sim.beta.cap[i,temp.cap.loc] * temp.covars[temp.cap.loc]))
        }else {
          temp.cap <- exp(sim.beta.zero.cap[i,1])
        }
        #Calculate Survival
        temp.surv <- temp.prod/(1+((temp.prod*sim.N[s+1,i,c])/temp.cap)) 
        #Update Pop Matrix
        sim.N[s+2,i,c] <- sim.N[s+1,i,c] * temp.surv  #update population matrix
      }#next f
      
      #Ocean
      s <- 1
      for(s in 1:n.o.stages) {

        #Calculate Maturation Schedule
        sim.prob.mat <- vector(length=n.oAges)
        for(a in 1:n.oAges) {
            sim.prob.mat[a] <- 1/(1+exp(sim.mat.param[i,1]*(5-a)+log(1/0.99-1))) # Probability of maturing at age
        }
        
        #Productivity
        temp.prod.loc <- prod.ref[s+n.fw.stages,which(prod.ref[s+n.fw.stages,]!=0)]
        if(length(temp.prod.loc)>0) {
          temp.prod <- 1/(1+exp(-1*sim.beta.zero.prod[i,s+n.fw.stages] - sum(sim.ind.prod[i,temp.prod.loc] * sim.beta.prod[i,temp.prod.loc] * temp.covars[temp.prod.loc])))
          # temp.prod <- 1/(1+exp(-1*sim.beta.zero.prod[i,1+n.fw.stages] - sum(sim.beta.prod[i,temp.prod.loc] * temp.covars[temp.prod.loc])))
        }else {
          temp.prod <- 1/(1+exp(-1*sim.beta.zero.prod[i,s+n.fw.stages]))
        }
        #Capacity
        temp.cap.loc <- cap.ref[s+n.fw.stages,which(cap.ref[s+n.fw.stages,]!=0)]
        if(length(temp.cap.loc)>0) {
          temp.cap <- exp(sim.beta.zero.cap[i,2] + sum(sim.beta.cap[i,temp.cap.loc] * temp.covars[temp.cap.loc]))
        }else {
          temp.cap <- exp(sim.beta.zero.cap[i,2])
        }
        #Calculate Survival
        temp.surv <- temp.prod/(1+((temp.prod*sim.N[s+n.fw.stages+1,i,c])/temp.cap)) 
        
        #Update Pop Matrix
        #Continuous Version
        sim.N[s+n.fw.stages+2,i,c] <- sim.N[s+n.fw.stages+1,i,c] * exp(-1*( (sim.sel_bycatch[i,s]*sim.fmort[i,1]) + (-1*log(temp.surv)) )) * (1-sim.prob.mat[s])
        
        #Returns
        #Continuous Version
        sim.ret[i,c] <- sim.ret[i,c] + sim.N[s+n.fw.stages+1,i,c] * exp(-1*( (sim.sel_bycatch[i,s]*sim.fmort[i,1]) + (-1*log(temp.surv)) )) *(sim.prob.mat[s]) * (1-sim.hr)
        
        sim.rec[i,c] <- sim.rec[i,c] + sim.N[s+n.fw.stages+1,i,c] * exp(-1*( (sim.sel_bycatch[i,s]*sim.fmort[i,1]) + (-1*log(temp.surv)) )) *(sim.prob.mat[s])
        
      }#next o
    }#next c - covar
    
    #===============================
    #BASE
    #Initialize Spawners
    base.N[1,i] <- sim.spawners
    #Initialize Eggs
    base.N[2,i] <- sim.eggs
      
    #Freshwater
    s <- 1
    for(s in 1:n.fw.stages) {
        #Update base model survival
        base.prod <- 1/(1+exp(-1*sim.beta.zero.prod[i,s]))
        base.cap <- exp(sim.beta.zero.cap[i,1])
        base.surv <- base.prod/(1+((base.prod*base.N[s+1,i])/base.cap)) 
        base.N[s+2,i] <- base.N[s+1,i] * base.surv
      }#next f
      
    #Ocean
    s <- 1
    for(s in 1:n.o.stages) {
        #Update base model survival
        base.prod <- 1/(1+exp(-1*sim.beta.zero.prod[i,s+n.fw.stages]))
        base.cap <- exp(sim.beta.zero.cap[i,2])
        base.surv <- base.prod/(1+((base.prod*base.N[s+n.fw.stages+1,i])/base.cap)) 
        #Update Popn Matrix
        #Continuous Version
        base.N[s+n.fw.stages+2,i] <- base.N[s+n.fw.stages+1,i] * exp(-1*((sim.sel_bycatch[i,s]*sim.fmort[i,1]) + (-1*log(base.surv)))) * (1-sim.prob.mat[s])
        #Update Recruitment
        #Continuous Version
        base.ret[i] <- base.ret[i] + base.N[s+n.fw.stages+1,i] * exp(-1*((sim.sel_bycatch[i,s]*sim.fmort[i,1]) + (-1*log(base.surv)))) * (sim.prob.mat[s]) * (1-sim.hr)
        
        base.rec[i] <- base.rec[i] + base.N[s+n.fw.stages+1,i] * exp(-1*((sim.sel_bycatch[i,s]*sim.fmort[i,1]) + (-1*log(base.surv)))) * (sim.prob.mat[s])
      }#next o
      
  }#next i - sample from posterior

  #Recruits-per-spawner
  sim.RpS <- sim.rec/sim.spawners
  base.RpS <- base.rec/sim.spawners
    
  #Egg to Rec Survival
  sim.E2R <- sim.ret/sim.eggs
  base.E2R <- base.ret/sim.eggs
  
  #Percent differences
  pct.RpS <- array(data=0, dim=c(n.samp, n.covars), dimnames=list(c(1:n.samp), covar.names))
  pct.E2R <- array(data=0, dim=c(n.samp, n.covars), dimnames=list(c(1:n.samp), covar.names))
  
  c <- 1
  for(c in 1:n.covars) {
    pct.RpS[,c] <- (sim.RpS[,c]-base.RpS)/base.RpS*100
    pct.E2R[,c] <- (sim.E2R[,c]-base.E2R)/base.E2R*100
  }#next c
  
#   caterplot(as.mcmc(pct.E2R), collapse=TRUE, reorder=FALSE, 
#             quantiles=list(outer=c(0.025,0.975),inner=c(0.25,0.75)))
#   mtext('Percent Difference In Survival', side=1, line=2.5)
#   abline(v=0, col=rgb(1,0,0,alpha=0.5), lwd=2, lty=1)
  
  #PLOT IT OUT 
  pct.E2R.list <- melt(pct.E2R)
  names(pct.E2R.list) <- c('sim','var.name','value')
  
  var.stage <- covar.stage[as.numeric(pct.E2R.list$var.name)]
  
  pct.E2R.list <- data.frame(pct.E2R.list,var.stage)
  pct.E2R.list$var.name <- factor(pct.E2R.list$var.name, levels=rev(prod.covar.names), ordered=TRUE)
  pct.E2R.list$var.stage <- factor(pct.E2R.list$var.stage,levels=stage.names, ordered=TRUE)
  
  #Plot it Out
  tmp <- ggplot(pct.E2R.list, aes(var.name, value, fill=var.stage))
  tmp <- tmp + ylab('Percent Difference in Egg-to-Recruit Survival (%)') + xlab('Covariate')
  if(pop=='chena') {
    tmp <- tmp + ggtitle('Chena River')	
  }else {
    tmp <- tmp + ggtitle('Salcha')	
  }
  tmp <- tmp + geom_hline(yintercept=0, col='red')
  tmp <- tmp + scale_y_continuous(limits = quantile(pct.E2R.list$value, c(0.0001, 0.999)))
  # tmp <- tmp + geom_boxplot(width=0.25, lwd=0.5, outlier.shape=NA)#width=0.1, lwd=0.1
  tmp <- tmp + geom_boxplot(lwd=0.5, outlier.shape=NA)
  tmp <- tmp + coord_flip()
  plot(tmp)
  tmp <- tmp + theme(legend.position='none')
  plot(tmp)
  tmp <- tmp + geom_boxplot(width=0.25, lwd=0.5, outlier.shape=NA)
  tmp <- tmp + geom_violin(alpha = 0.5, lwd=0.1, scale='width') 
  plot(tmp)
  
  #With Violins
  tmp <- ggplot(pct.E2R.list, aes(var.name, value, fill=var.stage))
  tmp <- tmp + ylab('Percent Difference in Egg-to-Recruit Survival (%)') + xlab('Covariate')
  if(pop=='chena') {
    tmp <- tmp + ggtitle('Chena River')	
  }else {
    tmp <- tmp + ggtitle('Salcha')	
  }
  tmp <- tmp + geom_hline(yintercept=0, col='red')
  tmp <- tmp + scale_y_continuous(limits = quantile(pct.E2R.list$value, c(0.0001, 0.999)))
  # tmp <- tmp + geom_boxplot(width=0.25, lwd=0.5, outlier.shape=NA)#width=0.1, lwd=0.1
  # tmp <- tmp + geom_boxplot(lwd=0.5, outlier.shape=NA)
  tmp <- tmp + geom_boxplot(width=0.25, lwd=0.5, outlier.shape=NA)
  tmp <- tmp + geom_violin(alpha = 0.5, lwd=0.1, scale='width') 
  tmp <- tmp + coord_flip()
  plot(tmp)
  tmp <- tmp + theme(legend.position='none')
  plot(tmp)
  
  #PLOT THE EXPLICIT DIFFERENCE IN RpS
  dim(sim.RpS)
  length(base.RpS)
  
#   temp.df <- data.frame(base.RpS,sim.RpS)
#   names(temp.df)[1] <- 'BASE'
  list.rps <- melt(sim.RpS)
  names(list.rps) <- c('sim','var.name','value')
  var.stage <- covar.stage[as.numeric(list.rps$var.name)]
  
  list.rps <- data.frame(list.rps,var.stage)
  
  #Add Base List
  list.base <- data.frame(c(1:n.samp),'BASE',base.RpS,'BASE')
  names(list.base) <- names(list.rps)
  
  #Glue them together
  list.rps <- rbind(list.base,list.rps)
  
  list.rps$var.name <- factor(list.rps$var.name, levels=rev(c('BASE',prod.covar.names)), ordered=TRUE)
  list.rps$var.stage <- factor(list.rps$var.stage,levels=c('BASE',stage.names), ordered=TRUE)
  
  #Plot it Out
  tmp <- ggplot(list.rps, aes(var.name, value, fill=var.stage))
  tmp <- tmp + ylab('Recruits-per-Spawner') + xlab('Covariate')
  if(pop=='chena') {
    tmp <- tmp + ggtitle('Chena River')	
  }else {
    tmp <- tmp + ggtitle('Salcha')	
  }
  #Change colors
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }
  tmp <- tmp + scale_fill_manual(values=c("#999999", gg_color_hue(3)))
  
  tmp <- tmp + geom_hline(yintercept=0, col='red')
  tmp <- tmp + scale_y_continuous(limits = quantile(list.rps$value, c(0.00001, 0.995)))
  # tmp <- tmp + geom_boxplot(width=0.25, lwd=0.5, outlier.shape=NA)#width=0.1, lwd=0.1
  tmp <- tmp + geom_boxplot(lwd=0.5, outlier.shape=NA)
  tmp <- tmp + coord_flip()
  plot(tmp)
  tmp <- tmp + theme(legend.position='none')
  plot(tmp)
  tmp <- tmp + geom_boxplot(width=0.25, lwd=0.5, outlier.shape=NA)
  tmp <- tmp + geom_violin(alpha = 0.5, lwd=0.1, scale='width') 
  plot(tmp)
  
  #With Violins
  tmp <- ggplot(list.rps, aes(var.name, value, fill=var.stage))
  tmp <- tmp + ylab('Recruits-per-Spawner') + xlab('Covariate')
  if(pop=='chena') {
    tmp <- tmp + ggtitle('Chena River')	
  }else {
    tmp <- tmp + ggtitle('Salcha')	
  }
  #Change colors
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }
  tmp <- tmp + scale_fill_manual(values=c("#999999", gg_color_hue(3)))
  
  tmp <- tmp + geom_hline(yintercept=0, col='red')
  tmp <- tmp + scale_y_continuous(limits = quantile(list.rps$value, c(0.00001, 0.995)))
  tmp <- tmp + geom_boxplot(width=0.25, lwd=0.5, outlier.shape=NA)
  tmp <- tmp + geom_violin(alpha = 0.5, lwd=0.1, scale='width') 
  tmp <- tmp + coord_flip()
  plot(tmp)
  tmp <- tmp + theme(legend.position='none')
  plot(tmp)
  
  
#   boxplot(pct.RpS, las=2, ylim=quantile(pct.RpS, probs=c(0.01,0.99)))
#   abline(h=0, col=rgb(1,0,0,alpha=0.5))
#   boxplot(pct.E2R, las=2, ylim=quantile(pct.RpS, probs=c(0.01,0.99)))
#   abline(h=0, col=rgb(1,0,0,alpha=0.5))
  
  
  
  
  #Histogram of dirichlet weight
  par(mfrow=c(1,1))
  ggs_density(ggs(out.mcmc), family='dirichWeight')
  plotPost(out$BUGSoutput$sims.list$dirichWeight, xlab='dirichWeight')
  plotPost(out$BUGSoutput$sims.list$pind, xlab='Prior Probability for Indicator Variables (pind)')
 
   plotPost(out$BUGSoutput$sims.list$sigma.prod, xlab='sigma.prod')
  # plotPost(out$BUGSoutput$sims.list$q.juvi, xlab='Catchability of Juvenile Index of Abundance (q.juvi)')
  
  #SLOPE OF PRODUCTIVITY AT SEA
  plotPost(out$BUGSoutput$sims.list$slp.beta.zero.prod, xlab='slp.beta.zero.prod')
   
  #BYCATCH STUFF
  plotPost(out$BUGSoutput$sims.list$sigma.bycatch, xlab='sigma.bycatch')
  plotPost(out$BUGSoutput$sims.list$bycatch_dirichWeight, xlab='bycatch_dirichWeight')
  
  #Heirarchical parameters
  par(mfrow=c(2,1), oma=c(0,0,0,0), mar=c(5,2,1,2))
  plotPost(out$BUGSoutput$sims.list$mu.Fmort, xlab='mu.Fmort')
  plotPost(out$BUGSoutput$sims.list$sd.Fmort, xlab='sd.Fmort')
  # ggs_running(ggs(out.mcmc), family='sd.Fmort')
  # plotPost(out$BUGSoutput$sims.list$sel_bycatch[,1], xlab='1-Ocean sel_bycatch')
  ggs_density(ggs(out.mcmc), family='sel_bycatch')
  
  par(mfrow=c(n.oAges,1), oma=c(0,0,0,0), mar=c(5,2,1,2))
  for(i in 1:n.oAges) {
    plotPost(out$BUGSoutput$sims.list$sel_bycatch[,i], xlab=paste(i,'-Ocean sel_bycatch',sep=''), xlim=c(0,1))
  }
  
  par(mfrow=c(3,3))
  for(y in 1:n.init.years) {
    plotPost(out$BUGSoutput$sims.list$est.init.spawners[,y], xlab=paste('Init Year: ',init.years[y],sep=''))
  }#next y
  
  #Plot Fishing Mortality Rate
  F_bycatch.years <- (min(years)+stage.offset[n.fw.stages+1]):(max(years)+stage.offset[n.fw.stages+n.o.stages])
  bycatch.fit.years <- 1991:2010 #Years to which 
  loc.fit.years <- which(F_bycatch.years %in% bycatch.fit.years)
  
  par(mfrow=c(2,1), oma=c(2,0,1,1), mar=c(1,4,0,0))
  boxplot(out$BUGSoutput$sims.list$F_bycatch, outline=FALSE, col='blue', ylab='Fmort Bycatch', xlab='', xaxt='n')
  HR_bycatch <- 1-exp(-1*out$BUGSoutput$sims.list$F_bycatch)
  boxplot(HR_bycatch, outline=FALSE, col='blue', ylab='Harvest Rate Bycatch', xlab='Year', xaxt='n')
  axis(side=1, at=c(1:length(F_bycatch.years)), F_bycatch.years, las=2)
  
  #===============================================================================
  #Plot Bycatch Fit
  plotPost(out$BUGSoutput$sims.list$bycatch_dirichWeight, xlab='bycatch_dirichWeight')
  
  #Bycatch Fit
  par(mfrow=c(1,1), oma=c(0,0,0,0), mar=c(4,4,1,1))
  boxplot(out$BUGSoutput$sims.list$bycatch.pred, col=rgb(1,0,0, alpha=0.5), ylab='BSAI Chinook Bycatch', 
          xlab='Year', outline=FALSE, xaxt='n')
  axis(1, at=c(1:length(bycatch.fit.years)), labels=bycatch.fit.years,las=2)
  points(out$BUGSoutput$sims.list$bycatch.obs[1,], pch=23, bg=rgb(0,0,1, alpha=0.5), cex=2)
  legend('topleft', legend=c('Observed','Predicted'), title='Bycatch', pch=c(23,22), 
         pt.bg=c(rgb(0,0,1, alpha=0.5), rgb(1,0,0, alpha=0.5)), pt.cex=2)
  
  #Posterior Predictive Bycatch Fit
  par(mfrow=c(1,1), oma=c(0,0,0,0), mar=c(4,4,1,1))
  boxplot(out$BUGSoutput$sims.list$bycatch.post.pred, col=rgb(1,0,0, alpha=0.5), ylab='BSAI Chinook Bycatch', 
          xlab='Year', outline=FALSE, xaxt='n')
  axis(1, at=c(1:length(bycatch.fit.years)), labels=bycatch.fit.years,las=2)
  points(out$BUGSoutput$sims.list$bycatch.obs[1,], pch=23, bg=rgb(0,0,1, alpha=0.5), cex=2)
  legend('topleft', legend=c('Observed','Posterior Predictive'), title='Bycatch', pch=c(23,22), 
         pt.bg=c(rgb(0,0,1, alpha=0.5), rgb(1,0,0, alpha=0.5)), pt.cex=2)
  
  #Plot posterior
  y.lim <- c(0, max(apply(out$BUGSoutput$sims.list$bycatch.pred, c(2), quantile, probs=0.975)))
  plot(x=bycatch.fit.years, y=out$BUGSoutput$sims.list$bycatch.obs[1,], type='l', lwd=2, ylim=y.lim,
       col=rgb(0,0,1, alpha=0.5), las=2, ylab='BSAI Chinook Bycatch', xlab='', main='Bycatch Fit')
  points(x=bycatch.fit.years, y=out$BUGSoutput$sims.list$bycatch.obs[1,], pch=23, bg=rgb(0,0,1, alpha=0.5), cex=2)
  axis(side=1, at=bycatch.fit.years, labels=FALSE)
  
  polygon(x=c(bycatch.fit.years,rev(bycatch.fit.years)), y=c(apply(out$BUGSoutput$sims.list$bycatch.pred,c(2), quantile, probs=0.975),
                                                     rev(apply(out$BUGSoutput$sims.list$bycatch.pred,c(2), quantile, probs=0.025))), 
          col=rgb(1,0,0,alpha=0.25), border=FALSE)
  polygon(x=c(bycatch.fit.years,rev(bycatch.fit.years)), y=c(apply(out$BUGSoutput$sims.list$bycatch.pred,c(2), quantile, probs=0.75),
                                                     rev(apply(out$BUGSoutput$sims.list$bycatch.pred,c(2), quantile, probs=0.25))), 
          col=rgb(1,0,0,alpha=0.25), border=FALSE)
  
  lines(x=bycatch.fit.years, y=apply(out$BUGSoutput$sims.list$bycatch.pred,c(2), median), col='red', lwd=1)
  
  legend('topleft', legend=c('Observed','Posterior'), pch=c(23,22), 
         pt.bg=c(rgb(0,0,1, alpha=0.5), rgb(1,0,0, alpha=0.5)), pt.cex=2)
  mtext(side=1, 'Year', line=3.5)
  
  #Plot posterior PREDICTIVE
  y.lim <- c(0, max(apply(out$BUGSoutput$sims.list$bycatch.post.pred, c(2), quantile, probs=0.975)))
  plot(x=bycatch.fit.years, y=out$BUGSoutput$sims.list$bycatch.obs[1,], type='l', lwd=2, ylim=y.lim,
       col=rgb(0,0,1, alpha=0.5), las=2, ylab='BSAI Chinook Bycatch', xlab='', main='Bycatch Fit')
  points(x=bycatch.fit.years, y=out$BUGSoutput$sims.list$bycatch.obs[1,], pch=23, bg=rgb(0,0,1, alpha=0.5), cex=2)
  axis(side=1, at=bycatch.fit.years, labels=FALSE)
  
  polygon(x=c(bycatch.fit.years,rev(bycatch.fit.years)), y=c(apply(out$BUGSoutput$sims.list$bycatch.post.pred,c(2), quantile, probs=0.975),
                                                             rev(apply(out$BUGSoutput$sims.list$bycatch.post.pred,c(2), quantile, probs=0.025))), 
          col=rgb(1,0,0,alpha=0.25), border=FALSE)
  polygon(x=c(bycatch.fit.years,rev(bycatch.fit.years)), y=c(apply(out$BUGSoutput$sims.list$bycatch.post.pred,c(2), quantile, probs=0.75),
                                                             rev(apply(out$BUGSoutput$sims.list$bycatch.post.pred,c(2), quantile, probs=0.25))), 
          col=rgb(1,0,0,alpha=0.25), border=FALSE)
  
  lines(x=bycatch.fit.years, y=apply(out$BUGSoutput$sims.list$bycatch.post.pred,c(2), median), col='red', lwd=1)
  
  legend('topleft', legend=c('Observed','Posterior Predictive'), pch=c(23,22), 
         pt.bg=c(rgb(0,0,1, alpha=0.5), rgb(1,0,0, alpha=0.5)), pt.cex=2)
  mtext(side=1, 'Year', line=3.5)
  
  #===============================================================================
  #Plot Bycatch Agecomp Fit
  cols <- rainbow(n=n.oAges)
  par(mfrow=c(2,1), oma=c(2,2,2,1), mar=c(3,4,1,1))
  ac.obs <- apply(out$BUGSoutput$sims.list$bycatch_agecomp.obs,c(2,3), mean)
  prop.obs <- ac.obs/apply(ac.obs, 1, sum)
  barplot(t(prop.obs), col=cols, ylab='Observed', las=2)
  
  ac.pred <- apply(out$BUGSoutput$sims.list$bycatch_agecomp.pred,c(2,3), median)
  prop.pred <- ac.pred/apply(ac.pred,1,sum)
  
  barplot(t(prop.pred), col=cols, ylab='Predicted', las=2, names.arg=bycatch.fit.years)
  
  mtext('Age Composition Proportion', side=2, outer=TRUE, font=2, cex=1.25, line=0.5)
  mtext('Year', side=1, outer=TRUE, font=2, cex=1.25, line=0.5)
  
  par(mfrow=c(1,1), mar=c(4,4,1,1), oma=c(0,0,0,0))
  
  lim <- c(min(prop.pred, prop.obs), max(prop.pred, prop.obs))
  
  plot(x=NULL, y=NULL, xlab='Observed Proportion', ylab='Predicted Proportion', main=pop, xlim=lim, ylim=lim)
  a <- 1
  for(a in 1:n.oAges) {
    points(x=prop.obs[,a], y=prop.pred[,a], pch=21, bg=cols[a])
  }#next a
  segments(x0=-1, y0=-1, x1=2, y1=2, lty=2, lwd=2)
  legend('topleft', legend=oAges, pch=21, pt.bg=cols, title='Ocean Age')
  
  
  
  #PLOT COMPARISON OF OBSERVED AND PREDICTED SPAWN SIGMA
  par(mfrow=c(2,1))
  #Observed
  x.lim <- c(0, max(spawnAbundOE,out$BUGSoutput$sims.list$sigma.spawn))
  hist(spawnAbundOE, xlim=x.lim, ylab='Observed', main='', col='blue')
  hist(out$BUGSoutput$sims.list$sigma.spawn, xlim=x.lim, ylab='Predicted', main='', col='red')
  mtext(side=3, text='Sigma Spawn', outer=TRUE, font=2)
  
  ##########
    # Calcualte correlations in covariates
  par(mfrow=c(1,1), oma=c(6,6,1,1))
    beta.prod.df <- data.frame(out$BUGSoutput$sims.list$beta.prod)
    names(beta.prod.df) <- covar.names
    corplot(cor(beta.prod.df))
  
  #Plot Specific Bycatch
#   par(mfrow=c(n.o.stages,1), oma=c(0,0,0,0), mar=c(2,4,0,1))
#   dim(out$BUGSoutput$sims.list$bycatch.obs)
#   dim(out$BUGSoutput$sims.list$bycatch.pred)
#   
#   s <- 1
#   for(s in 1:n.o.stages) {
#     xpos <- boxplot(out$BUGSoutput$sims.list$bycatch.pred[,,s], ylab=paste('Ocean',s))
#     points(out$BUGSoutput$sims.list$bycatch.obs[1,,s], pch=23, bg=rgb(1,0,0, alpha=0.5))
#   }#next s
  
  #Plot Comparisons
  # ggs_pairs(ggs(out.mcmc), family='beta.zero', lower=list(continuous='density'))
  # ggs_pairs(ggs(out.mcmc), family='sigma', lower=list(continuous='density'))
  # # ggs_pairs(ggs(out.mcmc), family='beta.prod', lower=list(continuous='density'))
  
  # ggs_compare_partial(ggs(out.mcmc), family='beta.prod')
  # ggs_running(ggs(out.mcmc), family='beta.prod')
  # ggs_running(ggs(out.mcmc), family='beta.zero')
  
  # #Plot traceplot for indicators to evaluate mixing
  # ggs_traceplot(ggs(out.mcmc, family='^ind.prod\\[[1235]\\]')) + theme_fivethirtyeight()
  # ggs_traceplot(ggs(out.mcmc, family='^ind.prod\\[[6789]\\]')) + theme_fivethirtyeight()
  
  
  
  #==========================
#   list.beta.prod <- melt(out$BUGSoutput$sims.list$beta.prod)
#   list.ind.prod <- melt(out$BUGSoutput$sims.list$ind.prod)
#   names(list.beta.prod) <- c('sim','var.num','value')
#   var.name <- prod.covar.names[list.beta.prod$var.num]
#   var.stage <- covar.stage[list.beta.prod$var.num]
#   
#   list.beta.prod <- data.frame(list.beta.prod, var.name, var.stage)
#   list.beta.prod$var.name <- factor(list.beta.prod$var.name, levels=rev(prod.covar.names), ordered=TRUE)
#   list.beta.prod$var.stage <- factor(list.beta.prod$var.stage,levels=stage.names, ordered=TRUE)
#   list.beta.prod <- data.frame(list.beta.prod, list.ind.prod$value)
#   names(list.beta.prod)
#   
#   
#   #Plot it Out
#   tmp <- ggplot(list.beta.prod, aes(var.name, value, fill=var.stage))
#   tmp <- tmp + ylab('Effect') + xlab('Covariate')
#   if(pop=='chena') {
#     tmp <- tmp + ggtitle('Chena River')	
#   }else {
#     tmp <- tmp + ggtitle('Salcha')	
#   }
#   tmp <- tmp + geom_hline(yintercept=0, col='red')
#   tmp <- tmp + scale_y_continuous(limits = quantile(list.beta.prod$value, c(0.001, 0.999)))
#   # tmp <- tmp + geom_boxplot(width=0.25, lwd=0.5, outlier.shape=NA)#width=0.1, lwd=0.1
#   tmp <- tmp + geom_boxplot(lwd=0.5, outlier.shape=NA)
#   tmp <- tmp + coord_flip()
#   plot(tmp)
#   tmp <- tmp + theme(legend.position='none')
#   plot(tmp)
#   tmp <- tmp + geom_boxplot(width=0.25, lwd=0.5, outlier.shape=NA)
#   tmp <- tmp + geom_violin(alpha = 0.5, lwd=0.1, scale='width') 
#   plot(tmp)
  
  
  dev.off()
  
  #================================================================================================
  ##### MANUSCRIPT FIGURES #####
  #Create file
#   manuscript.dir <- paste('Manuscripts/',date(), sep='')
#   file.create(manuscript.dir)
  
  
  ########## MANUSCRIPT: Covariates with inclusion probability ########## 
  #GGPLOT - Beta prod
#   pdf(paste('Manuscript/Covariate_',pop,'.pdf', sep=''), height=5, width=10)
#   
#   list.beta.prod <- melt(out$BUGSoutput$sims.list$beta.prod) 
#   names(list.beta.prod) <- c('sim','var.num','value')
#   var.name <- prod.covar.names[list.beta.prod$var.num]
#   var.stage <- covar.stage[list.beta.prod$var.num]
#   
#   list.beta.prod <- data.frame(list.beta.prod, var.name, var.stage)
#   list.beta.prod$var.name <- factor(list.beta.prod$var.name, levels=rev(prod.covar.names), ordered=TRUE)
#   list.beta.prod$var.stage <- factor(list.beta.prod$var.stage,levels=stage.names, ordered=TRUE)
#   
#   #Determine which samples to remove p(b|ind.prod==0)
#   list.ind.prod <- melt(out$BUGSoutput$sims.list$ind.prod)
#   loc.remove <- which(list.ind.prod$value==0)
#   list.beta.prod <- list.beta.prod[-loc.remove,]
#   
#   #With Violins
#   tmp <- ggplot(list.beta.prod, aes(var.name, value, fill=var.stage))
#   tmp <- tmp + ylab('Effect') + xlab('Covariate')
#   tmp <- tmp + ylab('Effect') + xlab('Covariate')
#   tmp <- tmp + geom_hline(yintercept=0, col='red')
#   tmp <- tmp + scale_y_continuous(limits = quantile(list.beta.prod$value, c(0.001, 0.999)))
#   tmp <- tmp + geom_violin(alpha = 0.5, lwd=0.1, scale='width') 
#   tmp <- tmp + geom_boxplot(width=0.25, lwd=0.5, outlier.shape=NA)
#   tmp <- tmp + coord_flip()
#   tmp <- tmp + theme(legend.position='none', 
#                      # axis.text = element_text(size = 20),
#                      panel.grid.major = element_line(colour = "white"),
#                      panel.background = element_rect(fill = "light gray"))
#   
#   # plot(tmp)
#   
#   #=======================================
#   #Inclusion Probability
#   list.ind.prod <- melt(out$BUGSoutput$sims.list$ind.prod)
#   names(list.ind.prod) <- c('sim','var.num','value')
#   var.name <- prod.covar.names[list.ind.prod$var.num]
#   var.stage <- covar.stage[list.ind.prod$var.num]
#   
#   list.ind.prod <- data.frame(list.ind.prod, var.name, var.stage)
#   list.ind.prod$var.name <- factor(list.ind.prod$var.name, levels=rev(prod.covar.names), ordered=TRUE)
#   list.ind.prod$var.stage <- factor(list.ind.prod$var.stage,levels=stage.names, ordered=TRUE)
#   #Add the pind median to list
#   med.pind <- median(out$BUGSoutput$sims.list$pind)
#   list.ind.prod <- data.frame(list.ind.prod,med.pind)
#   
#   #Plot it Out
#   tmp.2 <- ggplot(list.ind.prod, aes(var.name, fill=var.stage, alpha=0.75))
#   tmp.2 <- tmp.2 + stat_summary_bin(aes(y = value), fun.y = "mean", geom = "bar")
#   # tmp.2 <- tmp.2 + ylim(0,1.1)
#   tmp.2 <- tmp.2 + scale_y_continuous(expand=c(0,0), limits=c(0,1.05))
#   tmp.2 <- tmp.2 + coord_flip()
#   tmp.2 <- tmp.2 + ylab('Inclusion Probability') + xlab('Covariate')
#   tmp.2 <- tmp.2 + theme(legend.position='none', axis.text.y=element_blank(),
#                          axis.title.y=element_blank(),
#                          # axis.text = element_text(size = 20),
#                          panel.grid.major = element_line(colour = "white"),
#                          panel.background = element_rect(fill = "light gray"))#,
#   # plot(tmp.2)
# 
#   # grid.arrange(tmp, tmp.2, ncol=2, nrow=1, widths=c(4,2))
#   
#   plot_grid(tmp, tmp.2, ncol=2, rel_widths=c(3,1), align='h')
  
  
  ##########
#   #Calcualte correlations in covariates
#   corplot(cor(out$BUGSoutput$sims.list$beta.prod))
  
#   covars
#   covar.offset
#   
  covars.df <- data.frame(covars)
  names(covars.df) <- covar.names
  cor(covars.df)
  
  n.yrs.cv <- nrow(covars.df)
  
  new.covars.mtx <- matrix(nrow=n.yrs.cv+max(covar.offset), ncol=n.covars)
  
  i <- 1
  for(i in 1:n.covars) {
    # print(i)
    temp.loc <- max(covar.offset)-covar.offset[i]+1
    new.covars.mtx[(temp.loc:(temp.loc+n.yrs.cv-1)),i] <- covars.df[,i]
    
  }
  
  new.covars.df <- data.frame(new.covars.mtx)
  names(new.covars.df) <- covar.names
  
  cor(new.covars.df, use='complete.obs')
  
  # dev.off()

  ########## MANUSCRIPT: Plot All Model Fits Together ########## 
  # par(mfrow=c(5,2))
  
  
  #================================================================================================
  ########## Save the Entire Workspace ########## 
  remove(chena,chena.2,salcha,salcha.2)
  save.image(file=paste('Workspace/', pop, '_', model,'.RData', sep=''))
  
  return(out)
  
}# End wrapper function

start <- date()

print(date())
# out.chena <- wrapper_func(pop='chena', sims=sims, thins=thins, model=model, run.parallel=TRUE)
print(date())
out.salcha <-  wrapper_func(pop='salcha', sims=sims, thins=thins, model=model, run.parallel=TRUE)

end <- date()

#==========================
print(paste('sims:', sims, 'thins:', thins))
print(paste('START:',start))
print(paste('END:',end))


#Dirichlet trial
# pdf('Plots/Dirichlet Fun.pdf', heigh=8, width=6)
# par(mfrow=c(3,1))
# boxplot((rdirichlet(1e3, alpha=c(1,1,1,1))), main='alpha=1,1,1,1')
# boxplot((rdirichlet(1e3, alpha=c(10,10,10,10))), main='alpha=10,10,10,10')
# boxplot((rdirichlet(1e3, alpha=c(100,100,100,100))), main='alpha=100,100,100,100')
# 
# dev.off()


############## A little movie fun ###############
#PLOT FIT IN LOG SPACE
# par(mfrow=c(1,1), oma=c(0,0,2,0), mar=c(3,4,0,1))
# 
# 
# #Number of samples from posterior
# n.samp <- out$BUGSoutput$n.sims
# samp.seq <- round(seq(from=1, to=n.samp, length.out=500),0)
# n.samp.seq <- length(samp.seq)
# 
# #Color scheme
# colorRampAlpha <- function(..., n, alpha) {
#   colors <- colorRampPalette(...)(n)
#   paste(colors, sprintf("%x", ceiling(255*alpha)), sep="")
# }
# 
# cols <- colorRampAlpha(c('orange','red'), n=n.samp.seq, alpha=0.35)
# 
# saveGIF({
# i <- 1
# for(i in 1:n.samp.seq) {
#   # s <- samp.seq[i]
# 
#   #Adult
#   spawn.plot.yrs <- 1987:2010
#   
#   ln.spawn.pred <- log(out$BUGSoutput$sims.list$spawn.pred)
#   ln.spawn.obs <- log(out$BUGSoutput$sims.list$spawn.obs)
#   
#   y.lim <- c(min(apply(ln.spawn.pred,c(2), quantile, probs=0.025), 
#                  apply(ln.spawn.obs,c(2), mean)),
#              max(apply(ln.spawn.pred,c(2), quantile, probs=0.975), 
#                  apply(ln.spawn.obs,c(2), mean)))
#   
#   plot(x=spawn.plot.yrs, y=ln.spawn.obs[1,], type='p', pch=21, bg=rgb(0,0,1, alpha=1), cex=2,
#        xlab='', ylab='ln(Spawning Abundance)', ylim=y.lim)
#   axis(side=1, at=spawn.plot.yrs, labels=FALSE)
#   axis(1, labels=FALSE)
#   
#   lines(x=spawn.plot.yrs, y=ln.spawn.obs[1,], col=rgb(0,0,1, alpha=0.5), lwd=3)  
#   
#   #Model
#   j <- 1
#   for(j in 1:i){
#     s <- samp.seq[j]
#     lines(x=spawn.plot.yrs, y=ln.spawn.pred[s,], col=cols[s], lwd=2)
#   }#next j
# }
# }, interval = 0.1, movie.name = "Spawn Abund Converge.gif", ani.width = 600, ani.height = 600)
# 
# 
getwd()


# #Looking at posterior correlation
# out.chena$BUGSoutput$sims.list$mat.param
# 
# mtx.cor <- cor(cbind(out.chena$BUGSoutput$sims.list$mat.param,
#       out.chena$BUGSoutput$sims.list$beta.zero.prod,
#       out.chena$BUGSoutput$sims.list$beta.zero.cap[,1],
#       out.chena$BUGSoutput$sims.list$sel_bycatch
#       ))
# corplot(mtx.cor)
# 
# 
# mtx.cor <- cor(cbind(out.chena$BUGSoutput$sims.list$mu.Fmort,
#                out.chena$BUGSoutput$sims.list$sd.Fmort,
#                out.chena$BUGSoutput$sims.list$sigma.bycatch)) Fix sigma.bycatch?

# cor(cbind(out.chena$BUGSoutput$sims.list$beta.zero.prod,out.chena$BUGSoutput$sims.list$beta.zero.cap)) #Fix capacities?


