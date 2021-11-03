library(secr)
library(mcmcOutput)
library(nimble)

source("Nimble_Functions.r") # required for nimble model


n <- colSums(fire.capthist[,,]) ## Count up all deer individuals captured at each trap
n <- colSums(n) ## Sum both sessions together

M = 500 ## augmented superpopulation

J = nrow(traps(fire.capthist)) ## Number traps

A = 945 ## mask area of the fire study

Bayes_mask <- as.matrix(Fire.mask) ## Out of secr format

mask_habitat_variables <- read.csv("Fire_mask_habitat_vars.csv") ## landscape variables 

## Next chunk is manipulating mask to get rid of annoying NAs and make square
Bayes_mask <- Bayes_mask[!is.na(mask_habitat_variables$water.dist.clean),] ## Take out NA cells

mask_habitat_variables <- mask_habitat_variables[!is.na(mask_habitat_variables$water.dist.clean),] ## Ditto for variables

Bayes_mask[Bayes_mask[,1]< (-267900),]  <- NA ## Now make a square mask

mask_habitat_variables <- mask_habitat_variables[!is.na(Bayes_mask[,1]),] ## Take out NA variables

Bayes_mask <- Bayes_mask[!is.na(Bayes_mask[,1]),] ## Ditto for cells


## Back to defining the data
nCell = nrow(Bayes_mask)  ## number of rows to 

constants <- list(M=M,J=J,nCell=nCell,A=A) ## These values don't change across MCMC iterations or internal loops

data.pois <- list(n=n, K = rep(2, 144), X=traps(fire.capthist), grid=Bayes_mask, dist.water.pix=mask_habitat_variables$water.dist.clean) ## Data (changes per location)

code.pois <- nimbleCode({
  
  for(c in 1:nCell) { ## Loop through landscape cells
    log(mu[c]) <- A.0 + A.dist.water*dist.water.pix[c]  ## modeling density as function of intercept and distance to water
    pp[c] <- mu[c]/EN ## per cell relative density
  }
  
  EN<- sum(mu[1:nCell]) ## Total density
  psi <- EN/M ## Resulting probabilty of inclusion of a given individual from superpopulation in realized population
  
  for(i in 1:M) {
    w[i] ~ dbern(psi) ## Bernoulli trial for each possible individual in superpopulation to be included in population
    S[i] ~ dcat(pp[1:nCell]) ## categorical distribution for likelihood of individual activity center placement
    gx[i]<- grid[S[i],1] ## individual's x location
    gy[i]<- grid[S[i],2] ## individual's y location
    
    d2[i,1:J]<- (gx[i]-X[1:J,1])^2+(gy[i]-X[1:J,2])^2 ## Distance of given individual from traps
    prob[i,1:J]<- g0 * exp(-d2[i,1:J]/2/sigma^2) * w[i] ## observation model
    
  }
  
  for(j in 1:J) { ## Loop through traps for observation model here
    Ptrap[j]<- 1-prod(1-prob[1:M,j]) ## per trap likelihood of captures across individuals
  }
  n[1:J] ~ dpois_by_row(Ptrap[1:J],K[1:J]) ## Poisson model of individual activity centers and their capture
  
  g0 ~ dunif(0,1) ## prior on detection probability
  sigma ~ dnorm(175,75)    ## prior on space use
  A.0 ~ dnorm(0, sd=0.1)   ## prior on intercept
  A.dist.water ~ dnorm(0, sd=0.1) ## prior on effect of distance to water
  N <- sum(w[1:M]) ## Sum up individuals in realized population
  D <- N/A ## Density
  
})
#-------------------------------------------------------------------------

## Defining initial values (generally uninformed here)
S.init <- sample(1:nCell, M, replace = TRUE)
w.init <- rbinom(M,1,0.5)
A.0.init <- runif(1, -1, 1)
A.dist.water.init <- runif(1, -1, 1)
mu.init <- runif(nCell, -1, 1)

# Initial values
inits.pois <- list(g0=0.2,sigma=160,S=S.init, w=w.init,mu=mu.init,A.0=A.0.init, A.dist.water=A.dist.water.init) ## Inits for building the nimble model

## create the model object
deer_uscr_mod <- nimbleModel(code=code.pois, constants=constants, data=data.pois, inits=inits.pois, check = FALSE) ## Build nimble model

Rmcmc_deer <- compileNimble(deer_uscr_mod, showCompilerOutput = F) ## Compile it

ModSpec_deer <- configureMCMC(deer_uscr_mod, onlySlice = TRUE) # Slice sampling

## add a binary state sampler for each w node
ModSpec_deer$removeSamplers(c("w"), print = FALSE)
Nodes.pois <- deer_uscr_mod$expandNodeNames("w")
for(Node in Nodes.pois) ModSpec_deer$addSampler(target = Node, type = "binary", print=FALSE) ## Adding 

ModSpec_deer$resetMonitors() ## Don't want to track results of everything
ModSpec_deer$addMonitors(c("N","A.0","A.dist.water","g0","sigma","D", "N")) ## Key parameters to track

Cmcmc.deer <- buildMCMC(ModSpec_deer) ## Build the MCMC

#
Cmodel.deer <- compileNimble(Cmcmc.deer, project = deer_uscr_mod, resetFunctions = T)  ## Compile

ni<- 100 ## MCMC iterations
nb<- 10 ## Burn in
nc<- 1 ## Chains
nt<- 1 ## Thinning (none here)

#--------------------
# Run the sampler

inits.deer <- function(){list(g0=0.2,sigma=rnorm(1,175,75),S=S.init,w=w.init,mu=mu.init,A.0=A.0.init, A.dist.water=A.dist.water.init)} ## Written so that different inits each time a chain is run

samp.deer <- runMCMC(Cmodel.deer, niter = ni, nburnin = nb, nchains = nc, thin = nt, inits=inits.deer,  ## Run the chain!
                     samplesAsCodaMCMC = TRUE)

pars.deer <- mcmcOutput(samp.deer) ## for working with mcmcOutput package 
summary(pars.deer) ## Inspect parameters

diagPlot(pars.deer) ## Look at trace plots and parameter distributions. Perhaps introduce logit stabilization of g0
