## N-mixture models in R ####

## If you haven't use the AHMbook package, install it now
install.packages("AHMbook")

## Load in the librarys for this exercise
library(unmarked)
library(AHMbook)
library(dplyr)
library(AICcmodavg)
library(ggplot2)
library(raster)

## Load in the Finnmark dataset from the AHMbook package. These data were collected in Norway from 2005 to 2008. The are the counts of birds (by species) within 100 m of the point within a 15 minute count period
data("Finnmark")

## The data are stored in a list, so we are going to rearrange the data into the format we need. First we'll pull out the site data
sites<-Finnmark$sites

## Next well pull out the count data
counts<-(Finnmark$counts)

## Because there are multiple species in multiple years, it's stored in an array. For simplicity, we'll select one year of data for one species. In this case, it's the Lapland Bunting in 2006
dimnames(counts)
lb_06<-counts[,,2,8]

## We might also be interested in the time of day that the points were surveyed 
time<-Finnmark$timeOfDay

## We want to select the year of data associated with the count data
dimnames(time)
time_2006<-time[,,2]

## We can now make our unmarked dataframe. Similar to occupany, we specify our count data (y), our site covariates, and our observation covariates
umf<-unmarkedFramePCount(y=lb_06, siteCovs = sites, obsCovs = list(time=time_2006))

## First we are going to figure out the mixture most suited for the data. We will try Poisson vs. Negative Binomial
null_p<-pcount(~1 ~1, data=umf, mixture="P")                         
null_nb<-pcount(~1 ~1, data=umf, mixture="NB")

mix_aic<-fitList(null_p, null_nb)
modSel(mix_aic)

## The negative binomial model fit better, so we will use that moving forward and next test models that affect detection while holding abundance constant. For this type of model, we specify detection covariates first in the formula, followed by abundance covariates

## A model where time of day affects detection
p_time<-pcount(~time ~1, data=umf, mixture="NB")
summary(p_time)

## A model where the area of willows in a plot affects detection
p_area<-pcount(~area ~1, data=umf, mixture="NB")
summary(p_area)

## We can now rank the detection models
p_aic<-fitList(null_nb, p_time, p_area)
modSel(p_aic)

### We'll retain the time of day covariate and now run models with covariates on abundance

## Area of willows affects abundance
abun_area<-pcount(~time ~area, data=umf, mixture="NB")
summary(abun_area)

## Length of habitat edges affects abundance
abun_edge<-pcount(~time ~edge, data=umf, mixture="NB")
summary(abun_edge)

## Height of willows affects abundance
abun_height<-pcount(~time ~height, data=umf, mixture="NB")
summary(abun_height)

## Density of willows affects abundance
abun_dens<-pcount(~time ~density, data=umf, mixture="NB")
summary(abun_dens)

## Rank the abundance models
abun_aic<-fitList(p_time, abun_area, abun_edge, abun_height, abun_dens)
modSel(abun_aic)

## The model with willow height has the lowest AIC. We can run a goodness of fit test on this model, this time using the NMix.gof.test function. We'll only do 20 simulations but ideally you'd want 1000+
Nmix.gof.test(abun_height, nsim=20, report=T)

## Even with a small number of simulations, there's nothing obvious to suggest a poor fit

## We can now make some predictive plots starting with abundance
newdat=data.frame(height=seq(from=min(sites$height), to=max(sites$height), by=0.1))

Predicted_abun=predict(abun_height,type="state",newdata=newdat, append=T)

ggplot(Predicted_abun, aes(y=Predicted, x=height))+
  geom_line()+
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.5)+
  ylab("Predicted abundance")+
  xlab("Willow height (cm)")+
  theme_classic()

## Now a plot for detection
newdat2<-data.frame(time=seq(from=min(time_2006, na.rm = T), to=max(time_2006, na.rm=T), by=1))

Predicted_det<-predict(abun_height, type="det", newdata=newdat2, append=T)

## The times are in military time expresses as thousands, so we are going to convert the number to make it in a more understandable format
Predicted_det$time2<-Predicted_det$time/100

ggplot(Predicted_det, aes(y=Predicted, x=time2))+
  geom_line()+
  geom_ribbon(aes(ymin=lower, ymax=upper),alpha=0.5)+
  ylab("Detection probability")+
  xlab("Hour of day")+
  theme_classic()


## The last thing we might want to do is get a total abundance number across our entire study area, or map out how abundance varies spatially. 

## The first thing we need to do is figure out our effective sampling area. These point counts were conducted in 100 m radius circles. This means the estimated abundance is associated with an area of 31415.9 square m (area of a circle = pi*r^2) or 0.03 km2

## Import raster. Note - this is fake data which is really elevation from my PhD study area. The scale of numbers matched up closely enough to make it work for this example
height_ras<-raster("fake_height_data.tif")

## Make sure the name is the same as what was in the model
names(height_ras)<-"height"

plot(height_ras)

## We can check to see how large our raster cells are
res(height_ras)

## We'll round and say that they are 30 by 30. This is an issue because each cell is only 900 sq m compared to the 31415 sq m of the sample point area. One options is to resample the raster to the correct cell size

## We are first going to aggregate our raster. The factor that we are using is the length of a cell equal to the plot area divide by the current length of cell
height2<-aggregate(height_ras, fact=sqrt(31415)/30)
plot(height2)

##Let's check
res(height2)
sqrt(31415)

## Close enough for this exercise - the rounding is a bit off but we'll ignore that

## In order for the prediction to run, the raster needs to be saved as a raster stack
height_stack<-stack(height2)

## We can now run the predict function on the raster
height_predict<-predict(abun_height, newdata=height_stack, type="state")
plot(height_predict)

## We can now calculate the predicted total population size in our fake study area, as well as the lower and upper confidence intervals
cellStats(height_predict$Predicted, stat="sum") ## This is the total estimate
cellStats(height_predict$lower, stat="sum")
cellStats(height_predict$upper, stat="sum")
