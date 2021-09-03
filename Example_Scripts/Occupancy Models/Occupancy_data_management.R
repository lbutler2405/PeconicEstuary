#### Occupancy Models for dataset 1987-2021 ####

library(vegan)
library(tidyverse)
library(mvabund)
library(modeest)
library(ggplot2)
library(reshape2)
library(grid)
library(gridExtra)
#library(ggvegan)
library(ggrepel)
library(factoextra)
library(ggforce)

library(DiversityOccupancy)

## Function for *mode*
# Create the function.

# .csv file containing % abundance
master.sheet <- read.csv("species_master_spreadsheet.csv")
head(master.sheet)

# .csv file containing real abundances
raw.spp <- read.csv("all_spp_raw.csv")
head(raw.spp)

# merge the two dataframes on 'id' -> final number of samples = 9523
merge.df <- merge(master.sheet, raw.spp, by = "id")
head(merge.df)

View(merge.df)
merge.df <- merge.df[,-2:-42]

#write.csv(merge.df, "D:/Stony Brook/Peconic Estuary/PeconicEstuary/Example_Data/1987-2012/merge.df.csv")

stations <- sort(unique(merge.df$Station), decreasing = FALSE)
stations

years <- sort(unique(merge.df$Year, decreasing = FALSE))
years

week <- sort(unique(merge.df$Week, decreasing = FALSE))
week

month <- factor(x = unique(merge.df$Month), levels = unique(merge.df$Month))
month

season <- factor(x = unique(merge.df$Season), levels = unique(merge.df$Season))
season

merge.df.2 <- merge.df %>% arrange(Year) # arrange dataframe by year
merge.df.2 <- merge.df.2[,-17:-20]

merge.df.2 <- merge.df.2[, -1:-10]

dates <- merge.df.2[,1:6]
head(dates)
merge.df.2 <- merge.df.2[, -1:-6]

merge.df.2 <- merge.df.2[,-42] # drop grand total
merge.df.3 <- merge.df.2
merge.df.3[is.na(merge.df.3)] <- 0
head(merge.df.3)

spp.dat.yr <- NULL
spp.dat.yr <- cbind(dates$Year, dates$Station, dates$Date, dates$Season, merge.df.3)
names(spp.dat.yr)[1:4] <- c("Year", "Station", "Date", "Season")
spp.names <- colnames(master.sheet[,2:42])
names(spp.dat.yr)[5:45] <- spp.names
head(spp.dat.yr) ## Full dataframe with Year, Station and Species

names(merge.df)[21:61] <- spp.names
merge.df <- merge.df[,-62]

#### Change data to occupancy ####
st.yr <- spp.dat.yr[,1:4]
spp.occupancy <- spp.dat.yr[,-1:-4]
spp.occupancy <- spp.occupancy %>% mutate_if(is.numeric, ~1 * (. != 0))
head(spp.occupancy)

spp.occupancy <- cbind(st.yr, spp.occupancy)
head(spp.occupancy)

#############

#### Example using unmarked ####

#1. Prepare data
env.var <- merge.df[,1:20]
head(env.var)
spp.abund <- merge.df[,21:61]
head(spp.abund)

spp.occupancy <- spp.abund %>% mutate_if(is.numeric, ~1 * (. != 0))
head(spp.occupancy)

spp.occupancy[is.na(spp.occupancy)] <- 0
head(spp.occupancy)

# 2. Combine env. data with spp.data ####
data.all <- cbind(spp.occupancy, env.var)
head(data.all)

# 3. Filter by year and take 1 instance site from each site ####

x.for <- NULL
xx <- NULL

for (i in years) {
  xx <- data.all %>% 
    filter(Year == i)
  xx <- xx %>% arrange(Station)
  print(paste("Year ", i))
  
  for (j in stations){
    sort.s <- xx %>% arrange(Station)
    xy <- sort.s %>% 
      filter(Station == j)
    print(paste("Station ", j))
    xy <- data.frame(xy)
    xyx <- xy[1,]
    xyx$Year[is.na(xyx$Year)] <- i
    xyx$Station[is.na(xyx$Station)] <- j
    x.for <- rbind(x.for, xyx)
  }
}

View(x.for[,41:61])

x.for$Station <- rep(1:77, 25)
rownames(x.for) <- 1:nrow(x.for)

env.vars <- x.for[,42:ncol(x.for)]

temp.1987 <- x.for %>% filter(Year == 1987)
nrow(temp.1987)
View(temp.1987[41:61])

pivot.table <- pivot_longer(x.for, cols = 1:41)

#### Use anchovies as an example ####

anch.yrs <- data.frame(index = 1:77)
head(anch.yrs)

for (i in years){
  anch <- pivot.table %>% 
    filter(Year == i) %>%
    filter(name == "anchovy")
  col. <- data.frame(anch$value)
  anch.yrs <- cbind(anch.yrs, col.)
}

colnames(anch.yrs)[2:ncol(anch.yrs)] <- years

View(anch.yrs)

# Environmental variables for anchovies #

for (i in years){
  env <- x.for %>% 
    filter(Year == i) 
  env <- env[,42:ncol(env)]
  write.csv(env, paste('D:/Stony Brook/Peconic Estuary/PeconicEstuary/Example_Data/Occupancy/env_vars/first_quad_per_year/Year', i, '.csv'))
}

#### Sort out dataframes for occupancy model - ANCHOVIES ####

y.dat <- anch.yrs[,-1] 
head(y.dat) # has NAs

# Site Covariates do not change - assume it will be the start depth and end depth from 1987

sites <- x.for %>% filter (Year == 1987)
sites <- sites[,43:44]
head(sites)

siteCov <- sites

# Observation covariates for each year
# Create a dataframe for each covariate over the years

Surfs <- NULL # Surface Salinity

for (i in years){
  ss <- x.for %>% 
    filter(Year == i) 
  ss <- ss$SurfaceSalinity
  Surfs <- cbind(Surfs, ss)
}
colnames(Surfs) <- years
head(Surfs)

##
Botts <- NULL # Bottom Salinity

for (i in years){
  bs <- x.for %>% 
    filter(Year == i) 
  bs <- bs$BottomSalinity
  Botts <- cbind(Botts, bs)
}
colnames(Botts) <- years
head(Botts)

##
SurfDO <- NULL # Surface DO - not useable

for (i in years){
  sDO <- x.for %>% 
    filter(Year == i) 
  sDO <- sDO$SurfaceDO
  SurfDO <- cbind(SurfDO, sDO)
}
colnames(SurfDO) <- years
head(SurfDO)

##
BottDO <- NULL # Bottom DO

for (i in years){
  bDO <- x.for %>% 
    filter(Year == i) 
  bDO <- bDO$BottomDO
  BottDO <- cbind(BottDO, bDO)
}
colnames(BottDO) <- years
head(BottDO)

##
SurfTemp <- NULL # Surface Temp.

for (i in years){
  SurfT <- x.for %>% 
    filter(Year == i) 
  SurfT <- SurfT$SurfaceTemp
  SurfTemp <- cbind(SurfTemp, SurfT)
}
colnames(SurfTemp) <- years
head(SurfTemp)

##
BottTemp <- NULL # Bottom Temp.

for (i in years){
  BottT <- x.for %>% 
    filter(Year == i) 
  BottT <- BottT$BottomTemp
  BottTemp <- cbind(BottTemp, BottT)
}
colnames(BottTemp) <- years
head(BottTemp)

##
sech <- NULL # Bottom Temp.

for (i in years){
  sh <- x.for %>% 
    filter(Year == i) 
  sh <- sh$Secchi
  sech <- cbind(sech, sh)
}
colnames(sech) <- years
head(sech)

##
wk <- NULL # week

for (i in years){
  w <- x.for %>% 
    filter(Year == i) 
  w <- w$Week
  wk <- cbind(wk, w)
}
colnames(wk) <- years
head(wk)
wk <- data.frame(wk)
head(wk)
colnames(wk) <- years

wk[,1:ncol(wk)] <- lapply(wk[,1:ncol(wk)], factor)
summary(wk)

##
seas <- NULL # season

for (i in years){
  se <- x.for %>% 
    filter(Year == i) 
  se <- se$Season
  seas<- cbind(seas, se)
}
colnames(seas) <- years
head(seas)

obsCov <- list(SurfaceSal = Surfs,
               Bottomsal = Botts,
               BottomDO = BottDO,
               SurfaceT = SurfTemp,
               BottomT = BottTemp,
               Secchi = sech,
               Week = wk,
               Season = seas)

# Create an unmarked DataFrame
umf <- unmarkedFrameOccu(y = y.dat, siteCovs = siteCov, obsCovs = obsCov)
head(umf)

summary(umf)

# Standardise the data
umf@siteCovs$StartDepth <- scale(umf@siteCovs$StartDepth)
umf@siteCovs$EndDepth <- scale(umf@siteCovs$EndDepth)


umf@obsCovs$SurfaceSal <- scale(umf@obsCovs$SurfaceSal)
umf@obsCovs$Bottomsal <- scale(umf@obsCovs$Bottomsal)
umf@obsCovs$BottomDO <- scale(umf@obsCovs$BottomDO)
umf@obsCovs$SurfaceT <- scale(umf@obsCovs$SurfaceT)
umf@obsCovs$BottomT <- scale(umf@obsCovs$BottomT)
umf@obsCovs$Secchi <- scale(umf@obsCovs$Secchi)

# Fit the model assuming constant detection

fm <- occu(formula = ~1 ~1, data = umf)
fm

# assuming constant detection
fm1 <- occu(formula = ~ 1
            ~ StartDepth + EndDepth,
            data = umf)

fm1

# Only adding salinity and temperature
fm2 <- occu(formula = ~ SurfaceSal + Bottomsal + SurfaceT + BottomT
            ~ StartDepth + EndDepth,
            data = umf)

fm2


# Adding all observed covariates
fm3 <- occu(formula = ~ SurfaceSal + Bottomsal + SurfaceT + BottomT + BottomDO + Secchi + Season
            ~ StartDepth + EndDepth,
            data = umf)

fm3

## Model Selection ##

fit <- fitList('psi(.)p(.)' = fm1,
               'psi(SurfaceSal + Bottomsal + SurfaceT + BottomT)p(.)' = fm2,
               'psi(SurfaceSal + Bottomsal + SurfaceT + BottomT + BottomDO + Secchi + Season)' = fm3)

modSel(fit)

#### Proportion of Area Occupeid ###

# 1. Assuming perfect detection

siteValue <- apply(X = y.dat,
                   MARGIN = 1,
                   FUN = "max", na.rm = TRUE)

mean(siteValue)

# 2. Accounting for Imperfect Detection #

AICbest <- occu(formula = ~ SurfaceSal + Bottomsal + SurfaceT + BottomT + BottomDO + Secchi + Season
                ~ StartDepth + EndDepth,
                data = umf)

re <- ranef(AICbest)
EBUP <- bup(re, stat="mean")
CI <- confint(re, level=0.9)
rbind(PAO = c(Estimate = sum(EBUP), colSums(CI)) / 237)

#### Model Evaluation ####

fitstats <- function(fm3, 
                     method = "nonparboot") {
  
  observed <- getY(fm3@data)
  expected <- fitted(fm3)
  
  resids <- residuals(fm3,
                      method = "nonparboot")
  
  sse <- sum(resids^2,
             na.rm = TRUE)
  
  chisq <- sum((observed - expected)^2 / expected,
               na.rm = TRUE)
  
  freeTuke <- sum((sqrt(observed) - sqrt(expected))^2, 
                  na.rm = TRUE)
  
  out <- c(SSE = sse,
           Chisq = chisq,
           freemanTukey = freeTuke)
  
  return(out)
  
}

pb <- parboot(fm3,
              fitstats,
              nsim = 1000,
              report = TRUE,
              method = "nonparboot")

pb ## p-value greater than >> 0.05 is a good model.

par(mfrow = c(3,1))

plot(pb,
     main = "",
     xlab = c("SSE", "Chisq", "FT"))

#################

## Using the first two samples per quadrat per year ##

x.temp <- NULL
xx.temp <- NULL

for (i in years) {
  xx.temp <- data.all %>% 
    filter(Year == i)
  xx.temp <- xx.temp %>% arrange(Station)
  print(paste("Year ", i))
  
  for (j in stations){
    sort.temp <- xx.temp %>% arrange(Station)
    xy.temp <- sort.temp %>% 
      filter(Station == j)
    print(paste("Station ", j))
    xy.temp <- data.frame(xy.temp)
    xyx.temp <- xy.temp[1:2,]
    xyx.temp$Year[is.na(xyx.temp$Year)] <- i
    xyx.temp$Station[is.na(xyx.temp$Station)] <- j
    x.temp <- rbind(x.temp, xyx.temp)
  }
}

x.1987 <- x.temp %>% filter(Year==1987)
View(x.1987[,42:ncol(x.1987)])

x.1988 <- x.temp %>% filter(Year == 1988)
View(x.1988[,42:ncol(x.1988)])

x.1989 <- x.temp %>% filter(Year == 1989)
View(x.1989[,42:ncol(x.1989)])



#########

#### Multi-species occupancy model using DiversityOccupancy for the full dataset ####

library(DiversityOccupancy)

spp.occupancy
spp.occ <- subset(spp.occupancy, select= -c(Date, Season))

sites.occ <- read.csv("D:/Stony Brook/Peconic Estuary/PeconicEstuary/Example_Data/1987-2012/Site_Cov.csv")
head(sites.occ)
  
obs.occ <- read.csv("D:/Stony Brook/Peconic Estuary/PeconicEstuary/Example_Data/1987-2012/Obs_Cov.csv")
head(obs.occ)
