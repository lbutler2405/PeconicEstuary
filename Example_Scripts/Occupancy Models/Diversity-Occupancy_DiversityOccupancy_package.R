#### Abundance-Occupancy ####
library(DiversityOccupancy)
library(knitr)

data("IslandBirds")
View(IslandBirds)
kable(head(IslandBirds[1:8]))

# Site covariates
data("siteCov")
head(siteCov)

# Detection Covariates
names(Daily_Cov)
head(Daily_Cov[[1]])

#### Fiting models for abundance and predicting alpha diversity ####

## If *dredge = TRUE* only the best model will be selected

BirdDiversity.nodredge <- diversityoccu(pres = IslandBirds, sitecov = siteCov,
                               obscov = Daily_Cov,spp =  5, form = ~ Day + Wind + Time ~ Elev + Wetland + Upland, dredge = FALSE)

BirdDiversity.nodredge$models[[2]]

BirdDiversity <- diversityoccu(pres = IslandBirds, sitecov = siteCov, # build a diversity-occupancy model
                               obscov = Daily_Cov,spp =  5, form = ~ Day + Wind + Time ~ Elev + Wetland + Upland, dredge = TRUE)

BirdDiversity
BirdDiversity$models[[2]] ## look at the model

#### Model selection for alpha diversity modeling ####

glm.BirdDiverse <- model.diversity(BirdDiversity, method = "g", delta = 2, squared = TRUE)
glm.BirdDiverse$Table

# The responseplot.diver function takes a modeldiversity object and one of the variables used to predict the alpha diversity index, 
# and makes a plot showing the response of the diversity index against the selected variable. 

k.plot <- responseplot.diver(glm.BirdDiverse, Elev) ## ggplot
k.plot

#### Selecting conservation areas based on alpha diversity and abundance of species of conservation concern ####

library(raster)
data(Birdstack) # load the birdstack data
plot(Birdstack) # plot the birdstack data

newproj <- '+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0'
library(rgdal)

Birdstack2 <- stack(projectRaster(Birdstack, crs=newproj))
plot(Birdstack2)

Selected.area <- diversity.predict(model = BirdDiversity, diverse = glm.BirdDiverse, new.data = Birdstack2, quantile.nth = 0.65, species =
                                     c(F,T,T,F,F), kml = TRUE, name = "Selected_Bird_Area.kml")

plot(Selected.area$diversity.raster)
plot(Selected.area$species)
plot(Selected.area$priority.area)

#### Using the *unmarked* package ####

data <- read.csv(system.file("csv", "widewt.csv", package = "unmarked"))
head(data)

# year data set
y <- data[ ,2:4]

# site data set
siteCovs <-  data[ ,5:7]

# date of observation dataset
obsCovs <- list(date = data[ ,8:10],
                ivel = data[ ,11:13])

# convert to an *unmarked occupancy* object to use with the *occu* model
umf <- unmarkedFrameOccu(y = y, siteCovs = siteCovs, obsCovs = obsCovs)

summary(umf)

## Standardise the data by transforming to a logit scale ##
umf@siteCovs$elev <- scale(umf@siteCovs$elev)
umf@siteCovs$forest <- scale(umf@siteCovs$forest)
umf@siteCovs$length <- scale(umf@siteCovs$length)

umf@obsCovs$date <- scale(umf@obsCovs$date)
umf@obsCovs$ivel <- scale(umf@obsCovs$ivel)

# Fitting an occu model

# Formula type:  ~ detection formula ~ occupancy formula #


fm <- occu(formula = ~ 1 
           ~ 1,
           data = umf)
fm

# to convert back to original scale
backTransform(fm, type = "state")

## Add covariates to the model syntax 

fm1 <- occu(formula = ~ 1 
            ~ forest + elev + length,
            data = umf)
fm1

fm2 <- occu(formula = ~ date + ivel + forest # detection formula
            ~ forest + elev + length, # occupancy formula
            data = umf)
fm2


#### Select the best model from the three previously built ####

fit <- fitList('psi(.)p(.)' = fm,
               'psi(forest + elev + length)p(.)' = fm1,
               'psi(forest + elev + length)p(date + ivel + forest)' = fm2)

modSel(fit)


#### Accounting for imperfect detection - realistic ####

AICbest <- occu(formula = ~ forest + ivel
                ~ elev,
                data = umf)

re <- ranef(AICbest)
EBUP <- bup(re, stat="mean")
CI <- confint(re, level=0.9)
rbind(PAO = c(Estimate = sum(EBUP), colSums(CI)) / 237)

## Predict for a different site

occuPred <- predict(occuGG,
                    type = "state",
                    newdata = ggPred,
                    na.rm = TRUE,
                    inf.rm = TRUE)

# Produce a level plot for the area under consideration
levelplot(Predicted ~ ggPred$x + ggPred$y,
          data = occuPred,
          col.regions = rev(terrain.colors(100)),
          at = seq(0,1,length.out=101))

library(rstan)
library(blmeco)
library(sads)
