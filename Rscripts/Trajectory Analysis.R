#### Trajectory Analysis ####

install.packages('ecotraj')
install.packages("smacof")
install.packages("vegclust")
install.packages("RColorBrewer")


library(ecotraj)
library(smacof)
library(vegclust)

#Description of sites and surveys
sites = c(1,1,1,1,2,2,2,2,3,3,3,3)
surveys=c(1,2,3,4,1,2,3,4,1,2,3,4)

# Raw data table
xy<-matrix(0, nrow=12, ncol=2)
xy[2,2]<-1
xy[3,2]<-2
xy[4,2]<-3
xy[5:6,2] <- xy[1:2,2]
xy[7,2]<-1.5
xy[8,2]<-2.0
xy[5:6,1] <- 0.25
xy[7,1]<-0.5
xy[8,1]<-1.0
xy[9:10,1] <- xy[5:6,1]+0.25
xy[11,1] <- 1.0
xy[12,1] <-1.5
xy[9:10,2] <- xy[5:6,2]
xy[11:12,2]<-c(1.25,1.0)
raw <- cbind(sites,surveys,xy)

#Build a Distance matrix
D = dist(xy, upper = FALSE) # euclidean distance
D

#### Begin Trajectory Analsysis ####
# To begin our analysis of the three trajectories, we display them in an ordination space, 
# using function trajectoryPCoA.

oldpar <- par(mar=c(4,4,1,1))
trajectoryPCoA(D, sites, surveys, traj.colors = c("black","red", "blue"), lwd = 2,
               survey.labels = T)
legend("topleft", col=c("black","red", "blue"), 
       legend=c("Trajectory 1", "Trajectory 2", "Trajectory 3"), bty="n", lty=1, lwd = 2)

par(oldpar)

#### Trajectory lengths ####
trajectoryLengths(D, sites, surveys)

#### Trajectory Angles ####
# Calculates the angles between consecutive segments
trajectoryAngles(D, sites, surveys) # The larger the angle value, the more is trajectory changing in direction

# In Angles, rho is the mean resultant length of circular statistics
# Rho is used to assess the degree of homogeniety of angle values: 1 = all angle are the same

#### Overall directionality of a ecosystem trajectory ####
trajectoryDirectionality(D, sites, surveys) # accounts for angles and lengths of trajectories

#### The relative position of each point of a trajectory ####
trajectoryProjection(D, 1:4, 1:4)
trajectoryProjection(D, 5:8, 5:8)
trajectoryProjection(D, 7, 1:4) # projection of 3rd state of trajectory of site 2 (row 7) onto trajectory of site 1 (1:4)

trajectoryProjection(D, 9:12, 1:4) # projection of point of trajectory of site 3 (rows 9:12) onto trajectory of site 1 (rows 1:4)
# The curved pathway of site ‘3’ projects its fourth point to the same relative position as its second point

#### Trajectory Convergence ####
# Study their symmetric convergence or divergence
# Tests of convergence based on the trend analysis of the sequences of distances between points of the two trajectories

trajectoryConvergence(D, sites, surveys, symmetric = TRUE) # Symmetric convergence
# Values of the statistic (‘tau’) larger than 0 correspond to trajectories that are diverging, 
# Values lower than 0 correspond to trajectories that are converging

trajectoryConvergence(D, sites, surveys, symmetric = FALSE)

#### Distances between segments and trajectories ####
# Calculation of distances between directed segments

Ds <- segmentDistances(D, sites, surveys)$Dseg
Ds
# Distances between segments are affected by differences in both position, size and direction
# Check which combination has the maximum distance

#### Display distances between segments in two dimensions using mMDS ####

mMDS <- smacof::mds(Ds)
mMDS

xret <- mMDS$conf
xret
oldpar <- par(mar=c(4,4,1,1))

plot(xret, xlab="axis 1", ylab = "axis 2", asp = 1, pch = 21, bg = c(rep("black", 3), rep("red", 3), rep("blue", 3)),
     xlim = c(-1.5, 1), ylim = c(-1, 1.5))
text(xret, labels = rep(paste0("s", 1:3), 3), pos = 1)
legend("topleft", pt.bg = c("black", "red", "blue"), pch = 21, bty = "n", 
       legend = c("Trajectory 1", "Trajectory 2", "Trajectory 3"))

par(oldpar)

#### Changing measure of distance ####

trajectoryDistances(D, sites, surveys, distance.type = "Hausdorff") # Hausdorff: equal to the maximum distance between directed segments
trajectoryDistances(D, sites, surveys, distance.type = "DSPD") # Directed Segment Path Distance: an average of distances between segments
# DSPD is symmetrized distance
# Change symmetry by adding symmetrization = NULL

trajectoryDistances(D, sites, surveys, distance.type = "DSPD", symmetrization = NULL)

#### Structural dynamics in permanent plots ####
# E.g. analyze the dynamics of 8 permanent forest plots located on slopes of a valley in the New Zealand Alps

data("avoca")
avoca_strat # stratifiedvegdata object
head(avoca_strat)
temp <- data.frame(avoca_strat)

head(avoca_sites) # 8 sites
head(avoca_surveys) # 9 surveys

# Use function vegdiststruct from package vegclust to *calculate distances* between forest plot states in terms of structure and composition
avoca_D_man <- vegclust::vegdiststruct(avoca_strat, method = "manhattan", transform = function(x){log(x+1)})

#### Display trajectories for 'avoca_strat' ####

oldpar <- par(mar=c(4,4,1,1))
par(mfrow=c(1,1))
trajectoryPCoA(avoca_D_man, avoca_sites, avoca_surveys, traj.colors = RColorBrewer::brewer.pal(8, "Accent"),
               axes=c(1,2), length= 0.1, lwd = 2)
legend("topright", bty="n", legend = 1:8, col = RColorBrewer::brewer.pal(8, "Accent"), lwd = 2)

#### Using mMDS to plot trajectories for Forest Abundance ####
mMDS.forest <- smacof::mds(avoca_D_man)
mMDS

oldpar <- oldpar
trajectoryPlot(mMDS.forest$conf, avoca_sites, avoca_surveys, traj.colors = RColorBrewer::brewer.pal(8, "Accent"), 
               axes = c(1,2), length = 0.1, lwd = 2)
legend("topright", bty="n", legend = 1:8, col = RColorBrewer::brewer.pal(8,"Accent"), lwd=2)

par(oldpar)


#### Select a specific trajectory to assess ####
## plotTrajDistDiam function is a separate Rscript that needs to be amended based on the dataset?

oldpar <- par(mfrow=c(1,2))
trajectoryPCoA(avoca_D_man, avoca_sites, avoca_surveys, selection = 1, # Trajectory 1
               length = 0.1, lwd = 2, survey.labels = T)
plotTrajDiamDist(1)

oldpar <- par(mfrow=c(1,2))
trajectoryPCoA(avoca_D_man, avoca_sites, avoca_surveys, selection = 2, # Trajectory 2
               length = 0.1, lwd = 2, survey.labels = T)
plotTrajDiamDist(2)

oldpar <- par(mfrow=c(1,2))
trajectoryPCoA(avoca_D_man, avoca_sites, avoca_surveys, selection = 4, # Trajectory 4
               length = 0.1, lwd = 2, survey.labels = T)
# The turn at timepoints 8 and 9 indicates a recruitment of saplings
plotTrajDiamDist(4)

#### Trajectory lengths, angles and overall directionality using the Forest Data ####
# Lengths
trajectoryLengths(avoca_D_man, avoca_sites, avoca_surveys)

# Angles
avoca_angl <- trajectoryAngles(avoca_D_man, avoca_sites, avoca_surveys)
avoca_angl

# Directionality
avoca_dir <- trajectoryDirectionality(avoca_D_man, avoca_sites, avoca_surveys)
avoca_dir

# Assess the relationship between *trajectory directionality* and mean resultant vector using angular information
avoca_rho <- trajectoryAngles(avoca_D_man, avoca_sites, avoca_surveys, all=TRUE)$rho
oldpar <- par(mar=c(4,4,1,1))
plot(avoca_rho, avoca_dir, xlab = "rho(T)", ylab = "dir(T)", type = "n")
text(avoca_rho, avoca_dir, as.character(1:8))

par(oldpar)

#### Distances between trajectories ####

# Calculate the resemblance between forest plot trajectories

trajectoryDistances(avoca_D_man, avoca_sites, avoca_surveys)
# Low value in distances implies similarities
# looked rather close in position in the PCoA ordination of Ω with all trajectories, so probably 
# it is position, rather than shape which has influenced this low value

# Using mMDS to assess similarities in ordination space\
avoca_D_traj_man <- trajectoryDistances(avoca_D_man, avoca_sites, avoca_surveys, distance.type = "DSPD", verbose = FALSE)
mMDS.for <- smacof::mds(avoca_D_traj_man)
mMDS.for

x<- mMDS.for$conf[,1]
y<- mMDS.for$conf[,2]

oldpar <- par(mar=c(4,4,1,1))
plot(x,y, type = "p", asp=1, xlab = paste0("Axis 1"),
     ylab = paste0("Axis 2"), col = "black", bg = RColorBrewer::brewer.pal(8, "Accent"), pch = 21)
text(x, y, labels = 1:8, pos = 1)

temp2 <- centerTrajectories(avoca_D_man, avoca_sites, avoca_surveys)
temp2

temp2.dist <- trajectoryDistances(temp2, avoca_sites, avoca_surveys)
temp2.dist
# useful in cases if want to focus on spatio-temporal interaction while discarding spatial patterns that are constant in time

#### Transformation of dissimilarities in community trajectory analysis ####
# Effect of square root on a simple directional trajectory

# We begin by defining the species data of the trajectory itself. 
# The dataset consists of four rows (i.e. surveys) and four columns (species).

sites = rep(1,4)
surveys=1:4
spdata = rbind(c(35,30,20,15),
               c(50,25,15,10),
               c(65,20,10,5),
               c(80,15,5,0))
# Use function vegdist from package vegan to calculate the Bray-Curtis coefficient
D <- vegan::vegdist(spdata, "bray")
is.metric(D)
D # This dissimilarity matrix is a metric, so no need for transformation for CTA

trajectoryPCoA(D, sites, surveys) # Bray-Curtis dissimilarity responds linearly to the proposed sequence of community dynamics

# calculate the geometric properties of the trajectory
trajectoryLengths(D, sites, surveys)
trajectoryAngles(D, sites, surveys)
trajectoryDirectionality(D, sites, surveys)

#  take the square root of the dissimilarity values, 
# Necessary to achieve a metric (and Euclidean) space in a more complex dataset

sqrtD <- sqrt(D)
sqrtD

trajectoryPCoA(sqrtD, sites, surveys)
# The transformation increases all dissimilarity values (because the original values are smaller than 1), 
# but the increase is larger for smaller values, 
# so the ratio between large dissimilarities and small dissimilarities decreases

trajectoryLengths(sqrtD,sites,surveys)
trajectoryAngles(sqrtD,sites,surveys)
trajectoryAngles(sqrtD,sites,surveys, all=TRUE)
trajectoryDirectionality(sqrtD,sites,surveys)

#### Effect of different transformations on more complex trajectories ####
#  simulated dynamics to build another trajectory with more species:
# 1. 50 time steps
# 2. 50 individuals

Nsteps = 50
CC = 50
Nreplace <- CC*0.05 # Nreplace is number of individuals to be replaced each time step (5%)

# Define the initial community vector
x <- c(0, 1, 0, 67, 1, 3, 0, 2, 2, 2, 1, 6, 2, 0, 0, 2, 5, 1, 6, 0)

# vector with the probabilities of offspring for each species according
poffspring <- c(0, 0, 0.002, 0.661 ,0 ,0, 0.037, 0.281, 0, 0, 0, 0.008, 0, 0, 0.005, 0.003, 0, 0, 0, 0)

# simulate the dynamics by sequentially applying stochastic deaths and recruitment
m <- matrix(0, nrow=Nsteps+1, ncol = length(x))
m

m[1,] <- x
m

for(k in 1:Nsteps) {
  pdeath <- x/sum(x) # equal probability of dying
  deaths <- rmultinom(1, Nreplace, pdeath)
  x <- pmax(x - deaths, 0)
  offspring <- rmultinom(1, Nreplace, as.vector(poffspring))
  x <- x + offspring
  m[k+1,] <- x
}

# how frequently (with respect to the simulated step) a sample of the community is taken

Sj <- seq(1, Nsteps+1, by = 4) # every 4 steps
Sj

mj <- m[Sj, ]
mj # species are columns, sampling points are rows

surveys <- 1:length(Sj)
surveys 

sites <- rep(1, length(Sj))
sites

# calculate the Bray-Curtis dissimilarity
D <- vegan::vegdist(mj, )
D
is.metric(D) # Check if some triangles obey the triangle inequality (depending on the simulation)

pcoa <- trajectoryPCoA(D, sites, surveys, selection = 1, length = 0.1, axes = c(1,2))
pcoa
# the global transformation of principal coordinates analysis (PCoA) with negative eigenvalue correction is performed

pcoaD <- dist(pcoa$points) 
pcoaD # trajectory has some twists derived from stochasticity in death and recruitment

# Using Bray-Curtis Similarity
sqrtD <- sqrt(D)
pcoaSqrt <- trajectoryPCoA(sqrtD, sites, surveys, selection=1, length = 0.1, axes = c(1,2))
pcoaSqrt

#  transform dissimilarities using metric multidimensional scaling (mMDS)
res <- mds(D, ndim = length(Sj)-1, type = 'interval')
res

mmsd <- dist(res$conf)
mmsd
mmsd.plot <- trajectoryPlot(res$conf, sites, surveys, selection = 1, length = 0.1, axes = c(1,2))

## While the three plots look different, the differences are not striking (besides rotation issues)
par(mfrow=c(2,2))
trajectoryPCoA(D, sites, surveys, selection = 1, length = 0.1, axes = c(1,2), survey.labels = T)
trajectoryPCoA(sqrtD, sites, surveys, selection=1, length = 0.1, axes = c(1,2), survey.labels = T)
trajectoryPlot(res$conf, sites, surveys, selection = 1, length = 0.1, axes = c(1,2), survey.labels = T)

par(mfrow = c(1,1))

## compare the stress of the global solutions

stress0(D,pcoaSqrt$points, type="interval") # Sqrt has highest strss values
stress0(D,pcoa$points, type="interval")
stress0(D,res$conf, type="interval")

# calculate geometric properties we are not limited by ordination plots and we can take into account all dimensions
anglesD <- trajectoryAngles(D, sites, surveys) # local/reference point
angles.sqrtD <- trajectoryAngles(sqrtD, sites, surveys)
angles.pcoa.D <- trajectoryAngles(pcoaD, sites, surveys)
angles.mmds.D <- trajectoryAngles(mmsd, sites, surveys)

df <- as.data.frame(rbind(anglesD, angles.sqrtD, angles.pcoa.D, angles.mmds.D))
row.names(df) <- c("local", "global.sqrt", "global.pcoa", "global.mmds")
round(df,2)
df

# Overall Directionality for reference (D) and all 3 transformations (sqrtD, pcoaD, mmsd)
trajectoryDirectionality(D, sites, surveys)
trajectoryDirectionality(sqrtD, sites, surveys)
trajectoryDirectionality(pcoaD, sites, surveys)
trajectoryDirectionality(mmsd, sites, surveys)

# square root transformation distorts the angles and overall directionality of trajectories on the space defined 
# by the percentage difference (alias Bray-Curtis) dissimilarity index

#  the less harmful global transformation is provided by metric Multidimensional Scaling, 
# but the need to embed distances in an Euclidean space of all three transformations implies a stronger requirement than being a metric, 
# and results in distortion