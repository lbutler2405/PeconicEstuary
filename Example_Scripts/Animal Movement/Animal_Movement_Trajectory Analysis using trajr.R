#### Trajectory Analysis using trajr ####
install.packages("trajr")
library("devtools")
devtools::install_github("JimMcL/trajr")

library(trajr)

## Trajr works with trajectories, where a trajectory is a simplification of a real path travelled by an animal, 
# and is a set of 2-dimensional spatial coordinates with a third temporal dimension

# Within trajr, trajectories are represented by R objects with class Trajectory. 
# The TrajFromCoords function is used to create a Trajectory object from a set of x-y coordinates, with optional times.
coords <- data.frame(x = c(1, 1.5, 2, 2.5, 3, 4), 
                     y = c(0, 0, 1, 1, 2, 1), 
                     times = c(0, 1, 2, 3, 4, 5))
coords

#  it is first necessary to transform the positions to a suitable spatial projection such as UTM 
# (possibly by using spTransform from the rgdal package)

LatLong <- data.frame(X = c(56.85359, 56.85478), Y = c(-118.4109, -118.4035))
names(LatLong) <- c("X","Y")
LatLong
 
# Convert it to a sp object
coordinates(LatLong) <- ~ Y + X # Longitude first
LatLong

# Add a coordinate reference system
proj4string(LatLong) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")

# Project using spTransform
Utm <- spTransform(LatLong, CRS("+proj=utm +zone=11 ellps=WGS84"))
Utm

# Change to a dataframe
Utm.df <- data.frame(Utm)
Utm.df

## Transform from Coordinates to Trajectory
trj <- TrajFromCoords(coords)
trj
plot(trj)

#### Smoothing Trajectorie ####
trj <- TrajGenerate(200, random = TRUE, angularErrorSd = .25)
trj

# Plot original trajectory
plot(trj, lwd = 1, lty = 1)

# Create a smoothed trajectory, filter order 3, length 31
smoothed <- TrajSmoothSG(trj, p = 3, n = 31)
# Plot it in slightly transparent red
lines(smoothed, col = "#FF0000A0", lwd = 2)

legend("topright", c("Original", "Smoothed"), lwd = c(1, 2), lty = c(1, 1), col = c("black", "red"), inset = 0.01)

#### Resampling Trajectories ####
trj <- TrajGenerate(10, stepLength = 2)
trj

# Plot original trajectory with dots at trajectory coordinates
plot(trj, lwd = 2)
points(trj, draw.start.pt = FALSE, pch = 16, col = "black", cex = 1.2)

# Resample to step length 1
# The process of resampling a trajectory to a fixed step length is called rediscretization
resampled <- TrajRediscretize(trj, 1)

# Plot rediscretized trajectory in red
lines(resampled, col = "#FF0000A0", lwd = 2)
points(resampled, type = 'p', col = "#FF0000A0", pch = 16)
legend("topright", c("Original", "Rediscretized"), col = c("black", "red"), 
       lwd = 2, inset = c(0.01, 0.02))

#### Resampling in cases of uneven trajectory points ####
# Generate trajectory with a point every 2 hours and highly variable speed (which equates to step length)
trj <- TrajGenerate(10, stepLength = 1, fps = .5, timeUnits = "hours", linearErrorSd = .8)
trj

# Plot original trajectory with dots at trajectory coordinates
plot(trj, lwd = 2)
points(trj, draw.start.pt = FALSE, pch = 16, col = "black", cex = 1.2)

# Resample to 1 hourly steps
resampled <- TrajResampleTime(trj, 1)

# Plot rediscretized trajectory in red
lines(resampled, col = "#FF0000A0", lwd = 2)
points(resampled, type = 'p', col = "#FF0000A0", pch = 16)
legend("topright", c("Original", "Resampled"), col = c("black", "red"), 
       lwd = 2, inset = c(0.01, 0.02))

#### Trajectory Analaysis ####
#### Analysing Speed ####

# The functions *TrajVelocity* and *TrajAcceleration* estimate velocity and acceleration as vectors at each point along a trajectory

trj <- TrajGenerate()
trj

plot(trj, lwd = 2)
points(trj, draw.start.pt = FALSE, pch = 16, col = "black", cex = 1.2)

# Smooth before calculating derivatives
smoothed <- TrajSmoothSG(trj, 3, 101)
smoothed
plot(smoothed)

# Calculate speed and acceleration
derivs <- TrajDerivatives(smoothed)
derivs

# Plot change-in-speed and speed
plot(derivs$acceleration ~ derivs$accelerationTimes, type = 'l', col = 'red', 
     yaxt = 'n',
     xlab = 'Time (s)',
     ylab = expression(paste('Change in speed (', m/s^2, ')')))
axis(side = 2, col = "red")
lines(derivs$speed ~ derivs$speedTimes, col = 'blue')
axis(side = 4, col = "blue")
mtext('Speed (m/s)', side = 4, line = 3)
abline(h = 0, col = 'lightGrey')

#### Straightness Index ####

# Generate some trajectories for use in examples
n <- 100
# Random magnitude of angular errors
angularErrorSd <- runif(n, 0, 2)

# Generate some trajectories with varying angular errors
trjs <- lapply(1:n, function(i) TrajGenerate(500, stepLength = 2, 
                                             angularErrorSd = angularErrorSd[i]))
trjs

# Rediscretize each trajectory to a range of step sizes
stepSizes <- c(1, 2, 10)
reds <- lapply(stepSizes, function(ss) lapply(1:n, function(i) TrajRediscretize(trjs[[i]], ss)))
reds

# Calculate straightness (D/L) for all of the rediscretized trajectories
ds <- sapply(reds, function(rtrjs) sapply(1:n, function(i) TrajStraightness(rtrjs[[i]])))

# Calculate alternate straightness (r) for all of the rediscretized trajectories
rs <- sapply(reds, function(rtrjs) sapply(1:n, function(i) Mod(TrajMeanVectorOfTurningAngles(rtrjs[[i]]))))

# Plot both indices on the same graph
plot(rep(angularErrorSd, 3), rs,
     pch = 16, cex = .8,
     col = c(rep('red', n), rep('blue', n), rep('darkgreen', n)),
     xlab = expression(sigma[Delta]), ylab = "Straightness",
     ylim = range(c(ds, rs)))
points(rep(angularErrorSd, 3), ds,
       pch = 3, cex = .8,
       col = c(rep('red', n), rep('blue', n), rep('darkgreen', n)))

legend("bottomleft", c(expression(italic(r)), "D/L", paste("Step length", stepSizes)), 
       pch = c(16, 3, 16, 16), 
       col = c("black", "black", "red", "blue", "darkgreen"), inset = 0.01)
