## Spatial population dynamics model

library(secr)
library(nimble)
library(mcmcOutput)

setwd("C:/Users/Owner/Desktop/WildlifeDataAnalysis/WildlifeDataAnalysis/5_pop_Dynamics_Nimble") ## Wherever you downloaded the file

## Trap file
traps_fire <- read.table("Fire_Traps.txt")
traps_fire <- as.matrix(traps_fire[,2:3])

xlim <- c(min(traps_fire[, 1]) - 1000, max(traps_fire[, 1]) + 1000) ## 1 km buffer for state space
ylim <- c(min(traps_fire[, 2]) - 1000, max(traps_fire[, 2]) + 1000) ## Same in y dimension
ntraps <- nrow(traps_fire) ## Number traps

K <- 2 ## Two occassions per primary session
M <- 500 ## Unrealistically large augmented 'superpopulation' size

## Data in deer by trap format
Yarr <- read.csv("Deer_by_trap_format.csv")

Yarr <- array(c(as.matrix(Yarr[,1:2]), as.matrix(Yarr[,3:4])), dim = c(nrow(Yarr), 2, 2)) ## Get it into array

## OK, now on to the Nimble code (the model itself)
nimble_oSCR_code <- nimbleCode({
  
  ## Open population SCR component
  
  mean.phi ~ dunif(0, .5) ## Prior on survival
  
  for (i in 1:M){
    for (t in 1:(T-1)){
      phi[i,t] <- mean.phi ## Probability of survival
    } #t
  } #i
  
  for (t in 1:T){
    gamma[t] ~ dunif(0, 1) ## Prior on probability of being alive at time t
  } #t
  
  for (t in 1:T){ ## Time specific
    p0[t] ~ dunif(0,1) ## Prior on baseline probability of detection
    sigma[t] ~ dunif(0,500) ## Prior on within primary session movement parameter
  }
  for (i in 1:M){
    
    #  Activity centers
    s[i,1] ~ dunif(xlim[1],xlim[2])
    s[i,2] ~ dunif(ylim[1],ylim[2])
    
    ## Observation model
    for(t in 1:T){
      d2[i,1:ntraps,t] <- pow(pow(s[i,1]-traps_fire[1:ntraps,1],2) + pow(s[i,2]-traps_fire[1:ntraps,2],2),1) ## distances between traps and activity centers
      
      for(k in 1:K){ ## loop through the two occasions 
        
        mu[i,k,1:ntraps,t] <- p0[t]*exp( - d2[i,1:ntraps,t]/(2*sigma[t]^2))*z[i,t] ## probability of detection 
        cp[i,k,1:ntraps,t] <- mu[i,k,1:ntraps,t]/(1+ sum(mu[i,k,1:ntraps,t]))  ## capture probability
        cp[i,k,ntraps2,t] <- 1-sum(cp[i,k,1:ntraps,t])   # last cell = not captured
        Ycat[i,k,t] ~ dcat(cp[i,k,1:ntraps2,t]) ## Categorical distribution
        
      }  # k
    }  # t
    
    ## Demographic process model
    z[i,1] ~ dbern(gamma[1]) ## Alive or dead at start 
    
    for (t in 2:T){
      # State process
      q[i,t-1] <- 1-z[i,t-1]    # Availability for recruitment 
      available[i,t-1]<- prod(q[i,1:(t-1)])
      mu2[i,t] ~ dbern(phi[i,t-1] * z[i,t-1])  ## Survive or not?
      mu3[i,t] <-  mu2[i,t] + gamma[t] * available[i,t-1] ## Probability of surviving + probability of being alive * availability for recruitment (0 or 1)
      z[i,t] ~ dbern(mu3[i,t]) ## If survive or get recruited, alive!
      
    } #t
  } #i
  
  ## Track recruits, if we want
  for (i in 1:M){
    recruit[i,1] <- z[i,1] ## Initial recruits
    for (t in 2:T){
      recruit[i,t] <- (1-z[i,t-1]) * z[i,t] ## Subsequent sessions
    } #t
  } #i
  
  ## recruitment rate
  for (t in 2:T){
    rec_rate[t] <- sum(recruit[1:M,t])/N[t-1]
  }
  
  ## Derive population parameters
  for (t in 1:T){
    N[t] <- sum(z[1:M,t])        # Population size
    B[t] <- sum(recruit[1:M,t])  # Number of entries
    D[t] <- N[t]/area            # Density
  } #t
  
})

## Now need some initial values 

zst <- c(rep(1, M/2), rep(0, M/2)) ## Let's start with half alive and half dead
zst <- cbind(zst, zst)

## Propose start locations for each individual in superpopulation in each session
Sst <- cbind(runif(M, xlim[1], xlim[2]), runif(M, ylim[1], ylim[2]))
Sst <- array(Sst, dim = c(M, 2, 2))

## Bundle the intial values in a named list
inits <- list(z = zst, sigma = c(300, 300), gamma = runif(2, 0, 1), s = Sst[,,1], p0=c(0.2, 0.2))

data <- list(Ycat = Yarr,  traps_fire=traps_fire) ## Bundle data in named list

ntraps2<-ntraps+1  ## ntraps 2 is the no capture value for traps in the data
area<-(xlim[2]-xlim[1])*(ylim[2]-ylim[1])/1000000 ## Study area in km^2

constants<-list(T=2, K=K, M = dim(Yarr)[1], ntraps=ntraps, ## Bundle constants in named list
                ntraps2=ntraps2, ylim=ylim, xlim=xlim, area=area)

nimble_deer_fire_dynamics <- nimbleModel(code=nimble_oSCR_code, constants=constants, data=data, inits=inits, check = TRUE)

nimble_deer_fire_dynamics_comp <- compileNimble(nimble_deer_fire_dynamics)

nimble_deer_fire_dynamics_comp_MCMC_conf <- configureMCMC(nimble_deer_fire_dynamics_comp)
nimble_deer_fire_dynamics_comp_MCMC_conf$resetMonitors()
nimble_deer_fire_dynamics_comp_MCMC_conf$addMonitors(c("N","D","gamma","mean.phi",
                                                       "sigma","p0", "rec_rate"))

nimble_deer_fire_dynamics_comp_MCMC_confB <- buildMCMC(nimble_deer_fire_dynamics_comp_MCMC_conf)

nimble_deer_fire_dynamics_comp_MCMC_confB_C <- compileNimble(nimble_deer_fire_dynamics_comp_MCMC_confB, project = nimble_deer_fire_dynamics_comp,resetFunctions = FALSE)

nimble_deer_fire_dynamics_samps <- runMCMC(mcmc = nimble_deer_fire_dynamics_comp_MCMC_confB_C,niter=1000, nburnin=200, nchains = 2, samplesAsCodaMCMC = TRUE)

pars_nimble_deer_fire_dynamics <- mcmcOutput(nimble_deer_fire_dynamics_samps)
summary(pars_nimble_deer_fire_dynamics)

diagPlot(pars_nimble_deer_fire_dynamics)


