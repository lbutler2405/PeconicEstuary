#### Peconic Environmental Data ####
install.packages("tidyverse")
install.packages("mice")
install.packages("VIM")
install.packages("mvabund")
install.packages("VAST")
install.packages("vegan")
install.packages("rstan")
install.packages("blmeco")
install.packages("DiversityOccupancy")
install.packages("sf")
install.packages("dggridR")
install.packages("pdp")
install.packages("fitdistrplus")
install.packages("fields")
install.packages("devtools")
install_github('r-barnes/dggridR', vignette=TRUE) # problem installing
install.packages("sads")
install.packages("tweedie")

library(lubridate)
library(sf)
library(raster)
library(dggridR) ## Problem installing
library(pdp)
library(mgcv)
library(fitdistrplus)
library(viridis)
library(fields)
library(devtools)


install.packages("unmarked")

library(vegan)
library(tidyverse)
library(mvabund)

# Load the data

env.data <- read.csv("peconic_environmental_data.csv")
summary(env.data)

colSums(is.na(env.data)) # Number of missing values per column

env.data <- env.data[,-17] # Removes the 'Comments' columns

env.data <- subset(env.data, select = -c(SURF_DO100,BOT_DO100)) # removes the columns for DO100 surfance and bottom

## Initially, remove the ID and date but 
## leave the rest of the information in (including year and week)

env.data <- subset(env.data, select = -c(ID, DATE))


#### MICE Imputation
library(mice)
env.data.mt <- as.matrix(env.data)
md.pattern(env.data.mt)

library(VIM)
aggr_plot <- aggr(env.data.mt, col=c('navyblue','red'), numbers=TRUE, sortVars=TRUE, 
                  labels=names(env.data.mt), cex.axis=.7, gap=3, ylab=c("Histogram of missing data","Pattern"))

#### Imputing ####
temp.data <- mice(env.data.mt, m = 10, method = "pmm", maxit = 100, seed = 500)
completedData <- complete(temp.data,6)
View(completedData)


#### Filtering dataset
may.june <- env.data %>% 
  filter(Week >=19 & Week <= 28) %>%
  filter(Year >= 2009)
may.june
####

#### Ordination Analysis ####

#### NMDS ####
area_subdom.nmds <- metaMDS(area_subdom)
plot(area_subdom.nmds, type = "n", main = "NMDS subdominant species")
text(area_subdom.nmds, display = "species", cex = 0.7)
NMDS_sub <- data.frame(scores(area_subdom.nmds))
NMDS1_sub <- NMDS_sub$NMDS1
NMDS2_sub <- NMDS_sub$NMDS2

#### DCA ####
spp_abund_stand <- decostand(spp_abund, method = 'hellinger')
spp_abund_dca <- decorana(spp_abund_stand)
plot(spp_abund_dca, display = "sites", cex = 0.5)
plot(spp_abund_dca, display = "species", cex = 0.5)

#### CCA ####
abund_ca <- cca(spp_abund)
abund_ca_env <- cca(spp_abund ~ soil_pH + slope + Altitude + pct_water, data = ashtrees_env)
plot(abund_ca_env, display = "sites")
plot(abund_ca_env, display = "species")
plot(abund_ca_env, display = c("sites", "bp"), main = "CCA Ashtrees")
plot(abund_ca_env, display = c("species", "bp"))

#### Redundancy Analysis ####
#### RDA ####
env_factors <- data.frame(read.csv("Data/Biotic_data/ashtrees_env_variables_14-01-19.csv"))
env_factors[is.na(env_factors)] <- 0

area_dom_stand <- decostand(area_dom_matrix, method = "hellinger")
area_dom_cca <- cca(area_dom_stand)
area_dom_rda <- rda(area_dom_stand)
plot(area_dom_cca, display = "sites")
plot(area_dom_rda, display = "sites")

area_dom_rda.env <- rda(area_dom_stand ~ soil_pH + slope + Altitude + pct_water + length_10m + distance_10m + 
                          length_25m + distance_25m + length_35m + distance_35m + ditch_distance, 
                        data = env_factors, scale = FALSE) # Scaled = FALSE
plot(area_dom_rda.env, display = c("sites", "bp"), main = "RDA Dominant")
plot(area_dom_rda.env, display = c("species", "bp"))

#### Multivariate Abundnace Analysis ####
library(mvabund)

patchno_dom_mvabund <- mvabund(patchno_dom_spp)

# Produces a range of plots for 
#visualising multivariate abundance data and its relationship to environmental
#variables

plot.mvformula(log(patchno_dom_mvabund+1) ~ exp(env_var$soil_pH), main=" ",
               xlab="% Soil pH - Log Scale ", ylab="Abundance [log scale]",
               overall.main="Species Abundance vs soil pH",
               fg="grey", las=1, scale.lab="ss",t.lab="o", mfrow=c(4,3),log="x") 


boxplot(patchno_dom_mvabund, horizontal = TRUE, las = 2, main = "Number of Patches") # boxplot
meanvar.plot(patchno_dom_mvabund) # check mean variance

#### ManyGLM analysis ####
patchno_dom_glm_poisson <- manyglm(patchno_dom_mvabund ~ soil_pH + slope + pct_water 
                                   + Altitude, data = env_var, family = "poisson") # ?manyglm - negative binomial/poisson/binomial
drop1(patchno_dom_glm_poisson)#AICs for one-term deletions

## This could mean that one of our assumptions was wrong: either our mean-variance relationship was wrongly specified, 
## or our assumed relationship between our response and predictors was incorrect.
plot(patchno_dom_glm_poisson) # shoud not give a fan shape

resid_glm_area <- residuals(area_dom_glm_negbinomial)
plot(resid_glm_area)
summary(resid_glm_area)

summary(area_dom_glm_negbinomial, resamp="monte.carlo", test="wald", nBoot=300)

anova_patchno_dom_spp <- anova.manyglm(patchno_dom_glm_poisson, show.time = "all") # summary of species vs species (i.e. species used as variables)
capture.output(anova_patchno_dom_spp, file = "Results/DOMINANT/anova_glm_dom_spp.doc")

anova_ind_dom_spp <- anova(patchno_dom_glm_spp, p.uni="adjusted", show.time = "all") # try to get anova results for each species with each species in each quadrat
capture.output(anova_ind_dom_spp, file = "Results/DOMINANT/anova_ind_dom_spp.doc")

#### Generalised Additive Models: GAMs ####
library(mgcv)
library(tweedie)
library(statmod)

gam1 = gam(cover>0~s(elev),family=binomial,data=dat5)
summary(gam1)

#### MANYANY() ####
# Can use the *manyany* function 

#### MANYGAM ####
# ft=manyany("gam",abund,y~s(soil.dry),data=X,family="poisson")

## coral data example ##

data(tikus)
coral <- as.matrix(tikus$abund[1:20,])
sumSpp = apply(coral>0,2,sum)
coral <- coral[,sumSpp>6] ## cutting to just species with seven(!) or more presences to cut
## computation time. Maybe rerun with less (e.g. 4 or more presences) if curious and patient.
coralX <- tikus$x[1:20,]

ftTimeRep <- manyany("gam", coral, coral ~ time+rep, data=coralX,
                     family="poisson", var.power=1.5, composition=TRUE)

ftRep <- manyany("gam",coral, coral ~ rep, data=coralX,
                 family="poisson", var.power=1.5, composition=TRUE)

anova(ftRep,ftTimeRep,nBoot=9)

