#### Trajectory Analysis ####

library(ecotraj)
library(smacof)
library(vegclust)

# Use real data #
df <- read.csv("D:/Stony Brook/Peconic Estuary/Peconic_Project/Data/Modified_Data/Total_per_site.csv")
head(df)
df <- df %>% arrange(Station)

## For quadrat 1 over 25 years ##
df.2 <- df[1:25,]
View(df.2)

sites.2 <- df.2$Station
survey.2 <- df.2$Year

sites.2
survey.2

xy.2 <- df.2[,-1:-2]
head(xy.2)

#gat.dat.2 <- gather(df, key = 'Species', value = 'Abundance', 3:43)
#head(gat.dat.2)

xy.2.m <- as.matrix(xy.2)
View(xy.2.m)


#Build a Distance matrix - Example
D.2 = dist(xy.2, upper = FALSE) # euclidean distance
D.2
summary(D.2)

oldpar <- par(mar=c(4,4,1,1))
trajectoryPCoA(D.2, sites.2, survey.2, lwd = 2, traj.colors = "Blue",
               survey.labels = T)
legend("topleft", col=c("Blue"), 
       legend=c("Quadrat 1"), bty="n", lty=1, lwd = 2)

par(oldpar)

#Build a Distance matrix
D.man = vegdist(xy.2, method = "manhattan")
summary(D.man)

oldpar <- par(mar=c(4,4,1,1))
trajectoryPCoA(D.man, sites.2, survey.2, lwd = 2, traj.colors = "Blue",
               survey.labels = T)
legend("topleft", col=c("Blue"), 
       legend=c("Quadrat 1"), bty="n", lty=1, lwd = 2)
####

## Build for two quadrats ###
df.1.2 <- df %>% filter(Station %in% (1:2))
View(df.1.2)

sites.1.2 <- df.1.2$Station
survey.1.2 <- df.1.2$Year

sites.1.2
survey.1.2

xy.1.2 <- df.1.2[,-1:-2]
head(xy.1.2)
nrow(xy.1.2)

#gat.dat.2 <- gather(df, key = 'Species', value = 'Abundance', 3:43)
#head(gat.dat.2)

xy.1.2.m <- as.matrix(xy.1.2)
View(xy.1.2.m)


#Build a Distance matrix - Example
D.1.2 = dist(xy.1.2, upper = FALSE) # euclidean distance
D.1.2
summary(D.2)

oldpar <- par(mar=c(4,4,1,1))
trajectoryPCoA(D.1.2, sites.1.2, survey.1.2, lwd = 2, traj.colors = c("Blue", "green"),
               survey.labels = T)
legend("topleft", col=c("Blue" , "green"), 
       legend=c("Quadrat 1", "Quadrat 2"), bty="n", lty=1, lwd = 2)

par(oldpar)

#Build a Distance matrix
D.man.1.2 = vegdist(xy.1.2, method = "bray")
summary(D.man.1.2)

oldpar <- par(mar=c(4,4,1,1))
trajectoryPCoA(D.man.1.2, sites.1.2, survey.1.2, lwd = 1, traj.colors = c("Blue", "green"),
               survey.labels = T, cex = 0.7)
legend("bottomright", col=c("Blue", "green"), 
       legend=c("Quadrat 1", "Quadrat 2"), bty="n", lty=1, lwd = 1)

## Quadrat 1 and 2, side by side ##
par(mfrow = c(1,2))
trajectoryPCoA(D.man.1.2, sites.1.2, survey.1.2, lwd = 1, traj.colors = "Blue", selection = 1,
               survey.labels = T, length = 0.1, axes = c(1,2))
trajectoryPCoA(D.man.1.2, sites.1.2, survey.1.2, lwd = 1, traj.colors = "green", selection = 2,
               survey.labels = T, length = 0.1, axes = c(1,2))

## Using mMDS for Quadrat 1 and 2 ##

mMDS = smacof::mds(D.man.1.2)
mMDS

par(mfrow = c(1,1))
oldpar  <- par(mar=c(4,4,1,1))
trajectoryPlot(mMDS$conf,  sites.1.2, survey.1.2,
               traj.colors = c("Black", "grey"), 
               axes=c(1,2), length=0.1, lwd=2, survey.labels = T)
legend("bottomright", bty="n", legend = c("Quadrat 1", "Quadrat 2"), col = c("Black", "grey"), lwd=2)

temp <- cbind(mMDS$conf, sites.1.2, survey.1.2)
temp <- as.data.frame(temp)

#ggplot(temp, aes(D1, D2))+
#  geom_path(arrow = arrow())+
#  geom_point(data = temp[survey.1.2,])+
#  scale_y_reverse()+
#  scale_x_reverse()+
#  coord_fixed()

## For all quadrats in 1 year ## - does not make sense
q.1987 <- df %>% 
  filter(Year == "1987")
q.1987

sites.1987 <- q.1987$Year
survey.1987 <- q.1987$Station

q.1987 <- q.1987[,-1:-2]

D.1987 = vegdist(q.1987, method = "bray")
D.1987
summary(D.1987)

oldpar <- par(mar=c(4,4,1,1))
trajectoryPCoA(D.1987, sites.1987, survey.1987, lwd = 2, traj.colors = "Blue",
               survey.labels = T)
legend("topleft", col=c("Blue"), 
       legend=c("Quadrat 1"), bty="n", lty=1, lwd = 2)

#### Loop to generate a trajectory analysis for each quadrat over 25 years - using bray-curtis dissimlarity index ####
for (i in unique(df$Station)) {
  filt <- df %>% filter(Station == i)
  sites.x <- filt$Station
  survey.x <- filt$Year
  print(head(sites.x))
  print(head(survey.x))
  
  filt.x <- filt[,-1:-2]
  D.filt = vegdist(filt.x, method = "bray")
  
  
  oldpar <- par(mar=c(4,4,1,1))
  trajectoryPCoA(D.filt, sites.x, survey.x, lwd = 1, traj.colors = "grey",
                 survey.labels = T, cex = 0.5)
}

#### Loop to generate a trajectory analysis for each quadrat over 25 years - using metric multi-dimension scaling (mMDS) ####
par(mfrow = c(1,1))
for (i in unique(df$Station)) {
  filt <- df %>% filter(Station == i)
  sites.x <- filt$Station
  survey.x <- filt$Year
  print(head(sites.x))
  print(head(survey.x))
  
  filt.x <- filt[,-1:-2]
  print("....")
  
  dist.x <- vegdist(filt.x, method = "bray")
  D.filt.mds <- smacof::mds(dist.x)
  print("mMDS Performed")
  
  oldpar <- par(mar=c(4,4,1,1))
  pdf(paste("D:/Stony Brook/Peconic Estuary/Peconic_Project/Analysis/Plots/Trajectory_Analysis_1987-2012/Quadrat", i, ".pdf"), paper = "a4r",
      width = 10, height = 10)
  
  print(trajectoryPlot(D.filt.mds$conf, sites.x, survey.x, lwd = 1, traj.colors = "darkgrey",
                 survey.labels = T, length = 0.2))
  print(legend("bottomright", col="red", 
         legend=paste("S", i), bty="n", lty=1, lwd = 1, seg.len = 0))
  dev.off()
  
}

#### Trajectory Analysis for all quadrats BASED ON TOTAL PER SITE ####
head(df)
sites.f <- df$Station
survey.f <- df$Year
print(head(sites.f))
print(head(survey.f))

spp.f <- df[,-1:-2]
head(spp.f)
dist.f <- vegdist(spp.f, method = "bray")
mMDS.f <- smacof::mds(dist.f)
print("mMDS Performed")

par(cex = 0.8, cex.axis = 1.2, cex.lab = 1.2)
trajectoryPlot(mMDS.f$conf, sites.f, survey.f, lwd = 1, traj.colors = RColorBrewer::brewer.pal(5, "Accent"),
                     survey.labels = T, length = 0.2, cex = 2, selection = 1:5)
legend("bottomright", col=RColorBrewer::brewer.pal(5, "Accent"), 
             legend=paste("S", seq(1:5)), bty="n", lty=1, lwd = 1, seg.len = 0.5, pt.cex = 1, cex = 1.2)

#### Trajectory Analysis using average for each year, i.e. for the whole peconic ####

yearly <- NULL

for (i in unique(df$Year)) {
  annual <- df %>% filter(Year == i)
  annual <- annual[,-1:-2]
  annual <- colSums(annual)
  yearly <- rbind(yearly, annual)
}

rownames(yearly) <- unique(df$Year)
index <- rep(1, 77)
yearly <- cbind(years, index, yearly)
yearly <- as.data.frame(yearly)

# Trajectory Analysis
sites.y <- yearly$index
survey.y <- yearly$years
print(head(sites.y))
print(head(survey.y))

spp.y <- yearly[,-1:-2]
head(spp.y)
dist.y <- vegdist(spp.y, method = "bray")

mMDS.y <- smacof::mds(dist.y)
print("mMDS Performed")

par(cex = 0.8, cex.axis = 1.2, cex.lab = 1.2)
pdf("D:/Stony Brook/Peconic Estuary/Peconic_Project/Analysis/Plots/Trajectory_Analysis_1987-2012/Annual_combined/Stations_Combined_Bray.pdf", paper = "a4r",
    width = 10, height = 10)


trajectoryPCoA(dist.y, sites.y, survey.y, traj.colors = "grey", survey.labels = T, length = 0.2) # if using BRAY-CURTIS transformation

#trajectoryPlot(dist.y, sites.y, survey.y, lwd = 1, traj.colors = "grey", survey.labels = T, length = 0.2) # if using mMDS

dev.off()
