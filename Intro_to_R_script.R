###  Wildlife Population Analysis Intro to R

## Laura Gigliotti, Thomas Connor

# 10/12/2021

## First we need to check our working directory. This tells R where to access files and save files. 
# Your working directory should be the folder where you have stored your files/data

#Check your current working directory
getwd()

#If this isn't where your files are stored, change it to the correct filepath
#setwd("C:/Users/lgigliotti/Documents/Berkeley/Teaching/Intro")

## This code installs the packages - once they are installed you won't need to install them everytime you use R
#install.packages("tidyr")
#install.packages("dplyr")
#install.packages("ggplot2")

## This code loads the libraries - you will need to do this everytime you open R and want to use these packages
library(tidyr)
library(dplyr)
library(ggplot2)

# First we need to import our data. You should have a csv file called "state.csv" in the folder where your 
# R script is saved. We are going to name it "states" within R ("header=T" tells R that the first row in the 
# spreadsheet is the column names)

states<-read.csv("state.csv", header=T)

## Look at the structure of the data
str(states)

## We can also look at the first few rows of data using the "head" function
head(states)

## Similarly, we can look at the last few rows of data using the "tail" function
tail(states)

# This is an old dataset from the 70's that includes data on population sizes, income, 
# illiteracy rates, life expectancy rates, murder rates, high school graduation rates, 
# number of days below freezing, and area of US states 

## If we want to look at a specific column, we use the "$" to select the column
states$Name

## We can also look at values in specific cells by specifying the row and column number.
# We do this by defining the row number followed by the column number in a bracket. For example, 
# if we want to look at the 2nd row and 3rd column of our states dataframe we would use:
states[2,3]

## If we want to look at the entire 4th row:
states[4,]

## If we want to look at multiple rows or columns we can use a colon. For example, if we want 
# to look at columns 1 and 2:
states[,1:2]

## Finally, if we want to look at non-sequential rows or columns we can use c(). For example, 
# if we want to look at rows 2, 4, and 9:
states[c(2,4,9),]


## Let's look at some basic descriptive statistics of the data

## Average, minimum, and maximum life expectancy
mean(states$Life.Exp)
min(states$Life.Exp)
max(states$Life.Exp)

## Find the average, minimum, and maximum high school graduation rate
mean(states$HS.Grad)
min(states$HS.Grad)
max(states$HS.Grad)

##Sometimes you need to know what data minimum and maximum values are associated with
states[which.max(states$HS.Grad),]
states[which.min(states$HS.Grad),]

##You might be interested in knowing which rows meet a specific criteria (other than minimum or maximum).
# For example, we might want to know which states have a high school graduation rate greater than 60%
states[which(states$HS.Grad>60),]
states[which(states$Murder>15),]

##We can also use the "summary" function which will give us general summary statistics on all columns at once
summary(states$Income)

## We can add new columns using the $. Let's add a column for population density (population size divided by area)
states$Pop_density<-states$Population/states$Area

## Figure out the following:
 # 1) Which state has the highest murder rate?
states[which.max(states$Murder),] #Alabama
  
 # 2) Which state has the lowest illiteracy rate?
states[which.min(states$Illiteracy),] #Iowa
  
 # 3) Which states have illiteracy rates greater than 2%?
states[which(states$Illiteracy>2),] #AL, LA, MS, NM, SC, TX
  
## We can create some basic plots of our data. Maybe we are interested in the relationship between high school
# graduate rate and life expectancy. We will first make a plot using the standard R plotting function 
plot(Life.Exp~HS.Grad, data=states)

## This graph looks fine, but we need to add a title and change the axis labels. We might also be 
# interested in changing the color of the points

plot(Life.Exp~HS.Grad, data=states, main="Relationship between high school graduation rates and life expectancy in the US", 
     xlab="High school graduation rate", ylab="Life expectancy", col="blue")

## We can also make the same plot using a package called ggplot. This package lets you make graphs 
# with a lot more flexibility than base R.

ggplot(states, aes(x=HS.Grad, y=Life.Exp))+
  geom_point()+
  xlab("High school gradation rate")+
  ylab("Life expectancy")+
  ggtitle("Relationship between high school graduation rates and life expectancy in the US")

## We can also make different colors for different conditions. We will now separate points by region 
# in the United States. We will make the size of the points slightly bigger so that we can see them better.

ggplot(states, aes(x=HS.Grad, y=Life.Exp, color=Region))+ ##Adding the "color" identifier in this line adds different colors by regions
  geom_point(size=3)+ ## This line changes the size of the points
  xlab("High school gradation rate")+ ## This line adds an x-axis label
  ylab("Life expectancy")+ ## This line adds a y-axis label  
  ggtitle("Relationship between high school graduation rates and life expectancy in the US") ## This line adds a title

## Maybe we also want to label the points so that we know which state matches which point

ggplot(states, aes(x=HS.Grad, y=Life.Exp, color=Region))+ ##Adding the "color" identifier in this line adds different colors by regions
  geom_point(size=3)+ ## This line changes the size of the points
  xlab("High school gradation rate")+ ## This line adds an x-axis label
  ylab("Life expectancy")+ ## This line adds a y-axis label  
  ggtitle("Relationship between high school graduation rates and life expectancy in the US")+ ## This line adds a title
  geom_text(aes(label=Name), hjust=-0.1, vjust=-0.1)

## Finally, we might be interested in summarizing data based on categories. 

##For example, we might be interested in knowing how many states are in each region. We can use the "summarise" function to do this

state_count<-summarise(group_by(states, Region), count=n())

##We can also calculate different summary statistics such as the mean for each region. Our graph showed that there might be regional differences in life expectancy, so let's find the average life expectancy for each region

Region_lifeEx<-summarise(group_by(states, Region), Life_expectancy=mean(Life.Exp))

## We can now make a plot showing life expectancy by region. We will use the ggplot package once again.

ggplot(Region_lifeEx, aes(x= Region, y=Life_expectancy, fill=Region))+
  geom_bar(stat="identity")+
  ylab("Life expectancy")+
  ggtitle("Life expectancy in regions of the US")

## See if you can figure out how to make a figure showing how high school graduation rates differ by region

Region_Grad <- summarize(group_by(states,Region), HS.Grad=mean(HS.Grad))
ggplot(Region_Grad, aes(x = Region, y = HS.Grad, fill=Region))+
  geom_bar(stat="identity")+
  ylab("High School Graduation Rate")+
  ggtitle("High School Graduation Rate in Regions of the US")



## Let's transition to some more advanced R functions and models

## load needed packages
library(raster) ## Spatial package for working with raster data
library(dismo) ## Some nice functionality for species distribution modeling in this package
library(rgdal)
library(rgeos)

## let's read in some data

pres_data <- read.csv("Deer_scat_genotypes.csv") ## read em' in 

head(pres_data) ## OK, that is a lot of information (genetics, metadata)! But also Lat/Long in there

summary(pres_data$longitude) ## Looks like we have some NAs in the coordinates, don't want those!

pres_data <- pres_data[!is.na(pres_data$longitude),] ## Selecting rows that don't contain NA in the longitude column

pres_data_sp <- SpatialPointsDataFrame(coords = data.frame(pres_data$longitude, pres_data$latitude), data = pres_data, proj4string = CRS("+init=epsg:4326 +proj=longlat +ellps=WGS84
+datum=WGS84 +no_defs +towgs84=0,0,0"))  ## the "proj4string" is a projection that Thomas got off some of the Hopland spatial files
#Gotta specify projection and find it. If the file is already spatial, it has that info and you can just read in a shapefile

plot(pres_data_sp) ## Interesting ... what might be the cause of this pattern?


## Environmental Covariates
hopland_elevation <- raster("elevation.clean.tif") 
## You can also point to the file path directly when reading a file, in this case a raster depicting elevation values

plot(hopland_elevation)

plot(pres_data_sp, add=TRUE) ## They didn't show up, because they are not projected in UTMs!

pres_data_sp_utm <- spTransform(pres_data_sp, crs(hopland_elevation)) 
#crs argument matches the spatial projection of the first argument to the second
crs(hopland_elevation)#this shows you the projection information 

plot(pres_data_sp_utm, add=TRUE)


## Manipulating spatial data

## Let's load in a different raster dataset
RAP_data <- raster("RAP2018test.tif", band=5) ## This is a raster image with 6 bands, band 5 selects shrub cover. Edit file path as needed
#from the rangeland analysis platform that shows you veg cover estimates at 30m resolution

plot(RAP_data) ## Way too much spatial coverage!
points(pres_data_sp) ## Study area is a small subset of that

## OK, let's crop it to a km around the presence points
pres_data_buffer <- buffer(pres_data_sp, 1000) ## default unit meters, 1000 m is a 1km buffer
#defined buffer in lat/long, crop, and reproject into UTM because easier to reproject with smaller file

RAP_data_crop <- crop(RAP_data, pres_data_buffer) ## Crop to approx. study area

#Now we gotta reproject into UTMS
RAP_data_crop_utm <- projectRaster(from = RAP_data_crop, to = hopland_elevation) ## Convert to UTM, now that raster is smaller

RAP_data_crop_utm <- mask(RAP_data_crop_utm, hopland_elevation) ## Line up the layers nicely
#get rid of the areas we have no elevation data for (we're clipping to the elevation layer)
#crop only cuts to an extent
#mask just turns extra data into NAs
#generally, you should crop to get close and then mask to do something fine scale
plot(RAP_data_crop_utm)


deer_hab_predictors <- stack(hopland_elevation, RAP_data_crop_utm) ## Stack them

names(deer_hab_predictors) <- c("Elevation", "Shrub cover")
plot(deer_hab_predictors)


## Defining 'pseudoabsences' on landscape and preparing data

## Let's assume field techs had a 60 m search radius
survey_range <- buffer(pres_data_sp_utm, 60)
plot(survey_range)

elev_survey_mask <- mask(deer_hab_predictors$Elevation, survey_range)
plot(elev_survey_mask)
#for masking, first argument is the thing you're cropping, and survey range is what you crop to

pseudo_absenceses <- randomPoints(elev_survey_mask, n=1127, p=pres_data_sp_utm) ## Let's match number of presence points here (1127)

presence_covariates <- extract(deer_hab_predictors, pres_data_sp_utm) ## Extracts raster values at points, into matrix
presence_covariates <- as.data.frame(presence_covariates) ## Convert to data.frame 
presence_covariates$Presence <- 1 ## All these are presence locations

absence_covariates <- extract(deer_hab_predictors, pseudo_absenceses)
absence_covariates <- as.data.frame(absence_covariates)
absence_covariates$Presence <- 0 ## All these are absence locations

glm_data <- rbind(presence_covariates, absence_covariates) ## rowbind presence and absence data.frames together by common column names

## Fit a multivariate logistic regression
deer_presence_model <- glm(Presence ~ Elevation + Shrub.cover, data = glm_data, family = "binomial") ## Fit the model

summary(deer_presence_model) ## Summary of results

plot(deer_presence_model) ## Default is some model fit plots, these are really for linear models so not that useful


deer_presence_prediction <- predict(deer_hab_predictors,deer_presence_model, type = 'response') ## Predict to landscape

plot(deer_presence_prediction)#predicted probability of finding a scat in HOpland as function of elevation or shrub cover

## Not actually a 'population model' per se - what are we modeling?

## How to derive a 'population model' and predict relative abundance
## using this data?