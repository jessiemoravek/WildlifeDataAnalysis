## Multi-species Occupancy in R using the unmarked package

## Load packages
library(unmarked)
library(dplyr)
library(ggplot2)

## Load the first dataset. These data are camera trap detection/non-detection data from the AHMbook package (Rota et al. 2016). 
meso<-read.csv(file="occupancy_intro_data.csv", header=T)

## Load the second dataset
site_cov<-read.csv("occupancy_intro_sites.csv", header=T)

## We are going to look at how the occupancy of coyotes affects the occupancy of bobcats and vise versa. Let's filter out the data for each species.
coyote<-filter(meso, Species=="coyote")
bobcat<-filter(meso, Species=="bobcat")

## Similar to the simple occupancy model, we need to format our data into an unmarkedFrame. This time, our detection/non-detection data needs to be a list of data from each species in matrix format 

## We'll first save our detection data in the list/matrix format. Sp1 will be coyote and sp2 will be bobcat
y=list((as.matrix(coyote[,c(2:4)])),as.matrix(bobcat[,c(2:4)]))

## Here, we specify our detection list and our site covariate, similar to before
umf_multi<-unmarkedFrameOccuMulti(y=y,siteCovs=data.frame(People=site_cov$People_site, Dist=site_cov$Dist_5km, HDens=site_cov$HDens_5km, Trail=as.factor(site_cov$Trail)))

## Let's look at a summary of our unmarked frame
summary(umf_multi)
plot(umf_multi)                                  

## Running the model is similar to before, but the structure of the function is a bit different. The first argument are the detection covariates for species 1 followed by species 2. Note that the values need to be in quotes this time. The next argument is for the occupancy part of the model: covariates for species 1, covariates for species 2, and covariates for the intaction between the 2 species. In this first null model, we are going to set the interaction to 0 (e.g. species are occurring independently).

Null<-occuMulti(c("~1","~1"),c("~1","~1","~0"), data=umf_multi)
summary(Null)

## Similar to our simple occupancy model, we are first going to rank detection models to figure out the best detection model before testing models matching with occupancy hypotheses. For simplicity we are going to assume that detection covariates affect the two species similarly, but we could also have a larger suite of models where detection differs by species

## First we'll run a model where detection varies by # of people at a site
mod_det_people<-occuMulti(c("~People","~People"),c("~1","~1","~0"), data=umf_multi)

## And a model where detection varies depending if a camera is on a trail or not
mod_det_trail<-occuMulti(c("~Trail","~Trail"),c("~1","~1","~0"), data=umf_multi)

## Now we can rank our detection models
multi_det_aic<-fitList(Null, mod_det_people, mod_det_trail)
modSel(multi_det_aic)

## It looks like our trail detection model is the best, so well retain that. Now we get to the more complicated modeling of occurrence and co-occurrence.

## We'll first run a null model where we include the co-occurrence covariate (e.g. change the 0 from the last example to a 1)
mod_occ_null<-occuMulti(c("~Trail","~Trail"),c("~1","~1","~1"), data=umf_multi)

## Let's look at the results
summary(mod_occ_null)

## What we are interested in for these models is the "sp1:sp2" covariate. In this example, we see a significant positive coefficient - this means that coyotes are more likely to occupy sites where bobcats occur and vice versa. 

## There are lots of different combinations of covariates that we could include in models, but we are going to rank models to try to see what factors affect the cooccurrence of these species. Thus, we are only going to include 1 occupancy covariate per model and include it for both the occupancy and co-occurrence parts of the model

## First up, disturbance affects cooccurrence:
mod_occ_dist<-occuMulti(c("~Trail","~Trail"),c("~Dist","~Dist","~Dist"), data=umf_multi)
summary(mod_occ_dist)

## We can see an error message that the model didn't converge. Let's re-run the model and increase the number of iterations it's using to obtain estimates
mod_occ_dist<-occuMulti(c("~Trail","~Trail"),c("~Dist","~Dist","~Dist"), data=umf_multi, control=list(maxit=1000))
summary(mod_occ_dist)

## This time the model converged, but there are some large SE suggesting that we might need to take a closer look at this model later.

## In the next model, housing density affect occupancy/co-occurrence
mod_occ_house<-occuMulti(c("~Trail","~Trail"),c("~HDens","~HDens","~HDens"), data=umf_multi)
summary(mod_occ_house)

## Now we'll rank our occupancy models
multi_occ_aic<-fitList(mod_occ_null, mod_occ_dist, mod_occ_house)
modSel(multi_occ_aic)

## Housing density is our top model, so let's take a closer look at the relationships.
summary(mod_occ_house)

## We can see that we have covariates for the occupancy of coyotes and bobcats independently as well as the cooccurrence of both species. With these types of models we are typically most interested in conditional occupancy (what is the probably of occurrence of species A given that species B is present/absent). The best way to look at this is usually to make predictive plots. 

## Similar to before, we'll make a new dataframe based on the range of values of our covariate of interest
newdat=data.frame(HDens=seq(from=min(site_cov$HDens_5km), to=max(site_cov$HDens_5km), by=0.1))

# First we'll predict coyote occupancy when bobcat is present for the range of values
# in the newdata
coyote_bobcat <- predict(mod_occ_house, "state", species="sp1", cond="sp2", newdata=newdat)

## Let's add the original data values as well as a column indicating that bobcat is present (for some reason the predict function for this type of model doesn't recognize the "append" argument)
coyote_bobcat$HDens<-newdat$HDens
coyote_bobcat$Bobcat<-"Present"

## Now we'll do the same thing, but for sites where bobcat is absent. We do this by putting a minus sign in front of sp2
coyote_no_bobcat <- predict(mod_occ_house, "state", species="sp1", cond="-sp2", newdata=newdat)

## Let's add the original data values as well as a column indicating that bobcat is absent
coyote_no_bobcat$HDens<-newdat$HDens
coyote_no_bobcat$Bobcat<-"Absent"

## We can now combine our two predictions to make a figure
coyote_graph<-rbind(coyote_bobcat, coyote_no_bobcat)

ggplot(coyote_graph, aes(y=Predicted, x=HDens, color=Bobcat))+
  geom_line()+
  geom_ribbon(aes(ymax=upper, ymin=lower, fill=Bobcat), alpha=0.5)+
  theme_classic()+
  ylab("Occupancy probability")+
  xlab("Housing density within 5 km")


# NExt we'll predict bobcat occupancy when coyote is present for the range of values
# in the newdata
bobcat_coyote <- predict(mod_occ_house, "state", species="sp2", cond="sp1", newdata=newdat)

## Let's add the original data values as well as a column indicating that coyote is present
bobcat_coyote$HDens<-newdat$HDens
bobcat_coyote$Coyote<-"Present"

## Now we'll do the same thing, but for sites where coyote is absent. We do this by putting a minus sign in front of sp2
bobcat_no_coyote <- predict(mod_occ_house, "state", species="sp2", cond="-sp1", newdata=newdat)

## Let's add the original data values as well as a column indicating that coyote is absent
bobcat_no_coyote$HDens<-newdat$HDens
bobcat_no_coyote$Coyote<-"Absent"

## We can now combine our two predictions to make a figure
bobcat_graph<-rbind(bobcat_coyote, bobcat_no_coyote)

ggplot(bobcat_graph, aes(y=Predicted, x=HDens, color=Coyote))+
  geom_line()+
  geom_ribbon(aes(ymax=upper, ymin=lower, fill=Coyote), alpha=0.5)+
  theme_classic()+
  ylab("Occupancy probability")+
  xlab("Housing density within 5 km")
