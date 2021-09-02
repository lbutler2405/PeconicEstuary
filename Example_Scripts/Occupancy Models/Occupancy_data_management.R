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

write.csv(merge.df, "D:/Stony Brook/Peconic Estuary/PeconicEstuary/Example_Data/1987-2012/merge.df.csv")

stations <- sort(unique(merge.df$Station), decreasing = FALSE)
stations

years <- sort(unique(merge.df$Year, decreasing = FALSE))
years

week <- sort(unique(merge.df$Week, decreasing = FALSE))
week

month <- factor(x = unique(merge.df$Month), levels = unique(merge.df$Month))
month

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

#### Change data to occupancy ####
st.yr <- spp.dat.yr[,1:4]
spp.occupancy <- spp.dat.yr[,-1:-4]
spp.occupancy <- spp.occupancy %>% mutate_if(is.numeric, ~1 * (. != 0))
head(spp.occupancy)

spp.occupancy <- cbind(st.yr, spp.occupancy)
head(spp.occupancy)

#############

## Use anchovies as an example ##
anch <- NULL
anch <- spp.occupancy[,1:5]
head(anch)

# 1987
anch.1987 <- anch %>% filter(Year == 1987)
anch.1987 <- anch.1987 %>% arrange(Station)
anch.1987 <- anch.1987[,-1]
#date <- anch.1987$Date
#date <- as.character(date)
#date <- as.Date(date, format = "%d/%m/%Y")
#anch.1987$Date <- date
head(anch.1987)

View(anch.1987)
table(anch.1987$Season)

wide = anch.1987 %>% 
  spread(Station, anchovy)

View(wide)

#### Use total per year ###

tpy <- read.csv("D:/Stony Brook/Peconic Estuary/Peconic_Project/Data/Modified_Data/Total_per_site.csv")
head(tpy)
yr.st <- tpy[,1:2]
tpy <- tpy[,-1:-2]
tpy <- tpy %>% mutate_if(is.numeric, ~1 * (. != 0))
tpy <- cbind(yr.st, tpy)
head(tpy)

anch.tpy <- tpy[,1:3]
anch.tpy$anchovy[anch.tpy$anchovy > 0] <- 1 
head(anch.tpy)

wide = anch.tpy %>% 
  spread(Station, anchovy)

View(wide)

#### Loop for each species using total per year #### - problem occurs due to different sampling per site per year

tpy[,3:ncol(tpy)][,3:ncol(tpy) > 0] <- 1
head(tpy)

for (i in colnames(tpy[,3:ncol(tpy)])){
  spp <- tpy[,i]
  bind <- cbind(tpy[,1], tpy[,2], tpy[,i])
  bind <- as.data.frame(bind)
  colnames(bind) <- c("Year", "Station", paste(i))
  wide.i <- bind %>% 
    spread(Station, i)
  write.csv(wide.i, paste("D:/Stony Brook/Peconic Estuary/PeconicEstuary/Example_Data/Occupancy/", i, ".csv"))
  }

##############
df.temp <- merge.df[,-62]
env.temp <- merge.df[,1:20]
df.temp <- df.temp[,-1:-20]
colnames(df.temp) <- colnames(master.sheet[,2:42])
df.temp[is.na(df.temp)] <- 0
head(df.temp)

env.temp <- subset(env.temp, select = -c(StartLat, EndLat, StartLong, EndLong, Date))
df.temp <- cbind(env.temp, df.temp)
head(df.temp)

anch <- df.temp[,1:16]
anch <- anch[,-2:-10]
head(anch)

anch <- anch %>% arrange(Station)
head(anch)

anch.info <- NULL

# 1987 example # - anchovy

a.1987 <- anch %>% filter(Year == 1987)
head(a.1987)
a.1987 <- a.1987 %>% arrange(Station)

anch.info <- NULL
all.date <- NULL

for (i in unique(anch$Year)){
  lim <- anch %>% filter(Year == i)
  
  for (j in lim$anchovy){
    if (j > 0) {
      occ.1 <- 1
      anch.info <- rbind(anch.info, occ.1)
    } else {
      occ.0 <- 0
      anch.info <- rbind(anch.info, occ.0)
    }
    
    yr.st2 <- cbind(as.data.frame(i), as.data.frame(j))
    all.date <- rbind(all.date, yr.st2)
  }
}

anch.info
all.date

df.bind <- cbind(all.date, anch.info)
colnames(df.bind) <- c("Year", "Abundance", "Occupancy")
View(df.bind)

occupancy <- NULL

for (i in unique(df.bind$Year)) {
  yr <- df.bind %>% filter(Year == i)
  yr <- subset(yr, select = -c(Year, Abundance))
  #occupancy <- cbind(occupancy, yr)
}


###########

#### Ceate Observation list of dataframes #### - Use Year as columns, Site as rows and dataframe filled with week number?
# Example
spp.df <- read.csv("D:/Stony Brook/Peconic Estuary/PeconicEstuary/Example_Data/1987-2012/spp_df.csv")
head(spp.df)
spp.df <- spp.df %>% arrange(Year)
spp.df <- spp.df[,-45]
spp.df[is.na(spp.df)] <- 0

site.df <- read.csv("D:/Stony Brook/Peconic Estuary/PeconicEstuary/Example_Data/1987-2012/Site_Cov.csv")
site.df <- site.df %>% arrange(Year)
head(site.df)

obs.df <- read.csv("D:/Stony Brook/Peconic Estuary/PeconicEstuary/Example_Data/1987-2012/Obs_cov.csv")
obs.df <- obs.df %>% arrange(Year)
head(obs.df)

## 1987 example

y.1987 <- merge.df %>% filter(Year == 1987)
head(y.1987)
y.1987 <- y.1987[,-1:-10]
y.1987 <- y.1987[,-7:-ncol(y.1987)]
head(y.1987)
y.1987 <- y.1987 %>% arrange(Station)
y.1987 <- subset(y.1987, select = -c(Month, Date, Season))

observ <- NULL

for (i in unique(df$Year)) {
  xxl <- df.try %>% 
    filter(Year == i)
  xxxl <- xxxl %>% arrange(Station)
  xxl <- subset(xxl, select = -c(Month, Date, Season))
  print(paste("Year ", i))
}


