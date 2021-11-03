library(raster)
library(secr)

setwd("C:/Users/janel/Documents/Thesis/CamModels")

capthist_SMR <- read.capthist(captfile = "SMR_captfile.txt", 
                          trapfile = "SMR_trapfile_m_latlongfixed.txt", 
                          detector = "count", 
                          fmt = "trapID", 
                          trapcovnames = "attractant", 
                          markocc = c(rep(x = 0, times = 31)),
                          noccasions = 31
)

capthist_SMR <- addSightings(capthist_SMR, unmarked = "SMR_Tu.txt") ## y trap rows times x occasion columns

summary(capthist_SMR) ## Convenient summary

plot(capthist_SMR) # Plot

sightingPlot(capthist_SMR, mean = FALSE, dropunused = FALSE, scale = 50) ## Visualize counts of unmarked individuals

## Let's make the mask
mask_SMR <- make.mask(traps(capthist_SMR), buffer = 1000, spacing = 150) ## Same paramters as secr model

plot(mask_SMR) # plot to check shape

## 

buck_SMR <- secr.fit(capthist_SMR, mask = mask_SMR, detectfn = 0) ## secr fit function

buck_SMR ## Check results

# buck_SMR_2 <- secr.fit(capthist_SMR, mask = mask_SMR,start = buck_SMR, details = list(nsim = 1000)) ## Additional simulations to account for overdispersion, if needed
