#### Dynamic occupancy models in unmarked ####

install.packages("AICcmodavg")

M <- 250 # Number of sites
J <- 3 # num secondary sample periods
T <- 10 # num primary sample periods

psi <- rep(NA, T) # Occupancy probability
muZ <- z <- array(dim = c(M, T)) # Expected and realized occurrence
y <- array(NA, dim = c(M, J, T)) # Detection histories

set.seed(13973)

psi[1] <- 0.4 # Initial occupancy probability
p <- c(0.3,0.4,0.5,0.5,0.1,0.3,0.5,0.5,0.6,0.2)
phi <- runif(n=T-1, min=0.6, max=0.8) # Survival probability (1-epsilon)
gamma <- runif(n=T-1, min=0.1, max=0.2) # Colonization probability
# Generate latent states of occurrence
# First year
z[,1] <- rbinom(M, 1, psi[1]) # Initial occupancy state

# Later years
for(i in 1:M){ # Loop over sites
  for(k in 2:T){ # Loop over years
    muZ[k] <- z[i, k-1]*phi[k-1] + (1-z[i, k-1])*gamma[k-1]
    z[i,k] <- rbinom(1, 1, muZ[k])
  }
}

# Generate detection/non-detection data
for(i in 1:M){
  for(k in 1:T){
    prob <- z[i,k] * p[k]
    for(j in 1:J){
      y[i,j,k] <- rbinom(1, 1, prob)
    }
  }
}

# Compute annual population occupancy
for (k in 2:T){
  psi[k] <- psi[k-1]*phi[k-1] + (1-psi[k-1])*gamma[k-1]
}

plot(1:T, colMeans(z), type = "b", xlab = "Year",
     ylab = "Proportion of sites occupied",
     col = "black", xlim=c(0.5, 10.5), xaxp=c(1,10,9),
     ylim = c(0,0.6), lwd = 2, lty = 1,
     frame.plot = FALSE, las = 1, pch=16)
psi.app <- colMeans(apply(y, c(1,3), max))
lines(1:T, psi.app, type = "b", col = "blue", lty=3, lwd = 2)
legend(1, 0.6, c("truth", "observed"),
         col=c("black", "blue"), lty=c(1,3), pch=c(16,1))

# we reformat the detection/non-detection data from a 3-dimensional array (as generated) into a 2-dimensional matrix with M rows
yy <- matrix(y, M, J*T)
head(yy)

# create a matrix indicating the year each site was surveyed

year <- matrix(c('01','02','03','04','05','06','07','08','09','10'),
               nrow(yy), T, byrow=TRUE)
head(year)

## Sort the data to follow the *colext* requirement
simUMF <- unmarkedMultFrame(
  y = yy,
  yearlySiteCovs = list(year = year),
  numPrimary=T)

summary(simUMF)

# fit a few dynamic occupancy models
# 1.  fit a model with constant values for all parameters
m0 <- colext(psiformula= ~1, gammaformula = ~ 1, epsilonformula = ~ 1,
             pformula = ~ 1, data = simUMF, method="BFGS")
summary(m0)

# fit the dynamic occupancy model with full year-dependence in the parameters
# describing occupancy dynamics and also in detection

m1 <- colext(psiformula = ~1, # First-year occupancy
             gammaformula = ~ year-1, # Colonization
             epsilonformula = ~ year-1, # Extinction
             pformula = ~ year-1, # Detection
             data = simUMF)

m1

####  Manipulating results: prediction and plotting ####

## Below, we create data.frames called nd with
# each row representing a year. Then we request yearly estimates of the probability of
# extinction, colonization and detection, and compare them to “truth”, i.e., the values with
# which we simulated the data set



#### Real-type Dataset ####

library(ggplot2)
library(unmarked)
library(AICcmodavg)

# Load detection history (100 sites with 15 visits each)
detection_history <- read.csv("Dynamic_Occupancy_Intro/dynamic_detection_history.csv", row.names = "X") # presence-absence data
head(detection_history)

# Load covariate data
effort <- read.csv("Dynamic_Occupancy_Intro/dynamic_effort.csv",
                   # First variable ("X") has row.names but not data
                   row.names = "X") 
head(effort)

site_cov <- read.csv("Dynamic_Occupancy_Intro/dynamic_site_cov.csv",
                     # First variable ("X") has row.names but not data
                     row.names = "X")
head(site_cov)

#### Build an unmarkedmultFramOccu ####

frog_unmarkedMultFrame <- unmarkedMultFrame(y = as.matrix(detection_history), # 0's and 1's, one row per site, one column per survey
                                                          numPrimary = 5, #numPrimary is the number of primary surveys
                                                          obsCovs = list(effort=effort), 
                                                          siteCovs = site_cov[,2:4],
                                                          yearlySiteCovs = list(wetland_time = site_cov[,c(2, 5:8)]))
View(frog_unmarkedMultFrame)
summary(frog_unmarkedMultFrame)

#### Build the model ####

# initial occupancy (Ψ) is related to wetland amount (in the first year), temperature (temp), and precipitation (prec)
# detection probability (p) is related to search effort (scaled in our data to mean = 0, sd = 1)
# colonization probability (γ) is related to wetland amount (changes each year)
# extinction probability (ϵ) is related to wetland amount (changes each year)

dynamic_occ_m1 <- colext(
  psiformula = ~ wetland + temp + prec, # psi depends omn inital habitat, climate etc (initial environmental variables)
  gammaformula = ~ wetland_time, # colonization probability/depends on wetland change 
  epsilonformula = ~ wetland_time, # extinction probability/depends on wetland change
  pformula = ~ effort, # detection depends on survey effort
  data = frog_unmarkedMultFrame, 
  method = "BFGS") # quasi-Newton method 

dynamic_occ_m1

# Mackenzie-Bailey Goodness-of-Fit test
mb.boot <- AICcmodavg::mb.gof.test(dynamic_occ_m1, nsim = 20)
print(mb.boot, digit.vals = 4, digits.chisq = 4)

summary(dynamic_occ_m1) # compare p-values with the 'estimate' value. if p < 0.05 and estimate is +ve, then positively associated (and vice versa)

#### Smoothed Occupancy ####
# Get the smoothed *occupancy probabilities* for just site 25 (at year 1:5) 
# [2, , 25] = 2 for occupied, 
# " " for all 5 surveys, 25 for site 25;
# remove [] to print whole array
dynamic_occ_m1@smoothed[2, , 25]

# Mean smoothed occupancy probabilities can be accessed using:
smoothed(dynamic_occ_m1)

#### # Calculate SE for derived occupancy predictions using bootstrap ####
m1 <- nonparboot(dynamic_occ_m1, B = 20)
m1

# Predicted occupancy in each year (with SE)
# the "[2,]" calls the occupied estimates,
# "[1,]" for unoccupied estimates

predicted_occupancy <- data.frame(year = c(1:5), 
                                  smoothed_occ = smoothed(dynamic_occ_m1)[2,], 
                                  SE = m1@smoothed.mean.bsse[2,])
predicted_occupancy
plot(predicted_occupancy$smoothed_occ ~ predicted_occupancy$year)

plot <- ggplot(predicted_occupancy, aes(x=year, y=smoothed_occ)) + 
  geom_errorbar(aes(ymin=smoothed_occ-SE, ymax=smoothed_occ+SE), width=.1) +
  geom_line() +
  geom_point() +
  theme_bw()
plot + theme(axis.line = element_line(),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             panel.border = element_blank(),
             panel.background = element_blank()
             )
#### Comparing Models ####

# Fitting and comparing multiple occupancy models with the R package unmarked,
# Model-averaging predictions of occupancy*
# Model-averaging predicted relationships between occupancy and covariates

library(unmarked)
library(dplyr)
library(magrittr)
library(ggplot2)

# Load detection history (100 sites with 10 visits each)
detection_history <- read.csv("Dynamic_Occupancy_Intro/detection_history.csv", 
                              # First variable ("X") has row.names but not data
                              row.names = "X") 

# Examine data
head(detection_history)

# Load covariate data
effort <- read.csv("Dynamic_Occupancy_Intro/effort.csv",
                   # First variable ("X") has row.names but not data
                   row.names = "X") 
head(effort)

observers <- read.csv("Dynamic_Occupancy_Intro/observers.csv",
                      # First variable ("X") has row.names but not data
                      row.names = "X") 
head(observers)

site_cov <- read.csv("Dynamic_Occupancy_Intro/site_cov.csv",
                     # First variable ("X") has row.names but not data
                     row.names = "X")
head(site_cov)

# Add a random variable to site_cov (to demonstrate that not all variables should be included in the model)
# We'll assume this variable measures the amount of wetland in some buffer around the site
# Because it is a random number, it shouldn't affect occupancy in our example.
site_cov$wetland <- rnorm(n = nrow(site_cov), mean = 0, sd = 4)

head(site_cov)


#### Build an unmarkedFramOccu ####
sample.unmarkedFrame_cov <- unmarkedFrameOccu( 
  y = as.matrix(detection_history), # y is a matrix with observed detection history: 0's and 1's, one row per site, one column per survey
  obsCovs = list(effort = effort, observers = observers), # obsCovs = observation covariates in a list, each variable has site *rows x survey columns*
  siteCovs = site_cov) # siteCovs = dataframe with site rows x column variables

# S4 class for occupancy model data
summary(sample.unmarkedFrame_cov)

# Fit general model all variables
occu_p_full_psi_full <- occu(formula = ~effort + observers # detection formula first
                             ~forest + agri + wetland, # occupancy formula second,
                             data = sample.unmarkedFrame_cov)
occu_p_full_psi_full

dredge <- dredge(occu_p_full_psi_full)


#### Predictions using a model ####

library(AICcmodavg)
occu_modavg_psi_predict <- modavgPred(occu_p_full_psi_full, 
                                      # c.hat = 1, # to change variance inflation factor, default = 1) 
                                      parm.type = "psi", # psi = occupancy, can also be "det" for detection probability
                                      newdata = sample.unmarkedFrame_cov@siteCovs)[c("mod.avg.pred",
                                                                                     "lower.CL",
                                                                                     "upper.CL")]
