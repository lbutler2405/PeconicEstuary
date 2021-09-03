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

## Function for *mode*
# Create the function.
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

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

stations <- sort(unique(merge.df$Station), decreasing = FALSE)
stations

years <- sort(unique(merge.df$Year, decreasing = FALSE))
years

week <- sort(unique(merge.df$Week, decreasing = FALSE))
week

month <- factor(x = unique(merge.df$Month), levels = unique(merge.df$Month))
month

#### Grouping and Filtering ####
#### Number of times sample was taken per site per year ####
info <- NULL

for (i in unique(stations)){
  st.i <- merge.df %>% 
    filter(Station == i)
  table. <- table(st.i$Year)
  info <- dplyr::bind_rows(info, table.)
}

View(info)
info.df <- as.data.frame(info)
info.df[is.na(info.df)] <- 0

# Find mode for each year (columns)

mode.inf <- NULL # mode for each year

for (i in 1:ncol(info.df)){
  mode.i <- getmode(info.df[, i])
  mode.inf <- rbind(mode.inf, mode.i)
}

mode.inf <- t(mode.inf)
colnames(mode.inf) <- colnames(info)
mode.inf

info <- rbind(info, mode.inf)
row.names(info)[78] <- 'Mode'

min.yr <- apply(info.df, 2, min)
max.yr <- apply(info.df, 2, max)
mean.yr <- apply(info.df, 2, mean)
median.yr <- apply(info.df, 2, median)

calc.yr <- rbind(min.yr, max.yr, mean.yr, median.yr)

info <- rbind(info, calc.yr)
row.names(info)[79:82] <- c("Min", "Max", "Mean", "Median")
head(info)

## Find median, min, max and mean for each row

min.r <- data.frame(apply(info.df, 1, min))
max.r <- data.frame(apply(info.df, 1, max))
mean.r <- data.frame(apply(info.df, 1, mean))
median.r <- data.frame(apply(info.df, 1, median))

calc.rr <- cbind(min.r, max.r, mean.r, median.r)
colnames(calc.rr) <- c("Min", "Max", "Mean", "Median")
calc.rr[79:82, ] <- NA


info <- cbind(info, calc.rr)
#write.csv(info, "info_spp_yr.csv")

#### Number of times sample was taken per site ####
info.s <- NULL

st.y <- merge.df %>% 
  filter(Year == 1989)

head(st.y)
table. <- table(st.y$Station)
table.

for (i in unique(years)){
  st.y <- merge.df %>% 
    filter(Year == i)
  table. <- table(st.y$Station)
  info.s <- dplyr::bind_rows(info.s, table.)
}
rownames(info.s) <- years
View(info.s)

info.df <- as.data.frame(info)
info.df[is.na(info.df)] <- 0

##### Barplots - frequency of each site for each year ####

### Sample ####

index <- 1:nrow(info.df)
info.m <- cbind(index, info.df)
info.m <- info.m[-78:-82,-27:-30]


gat.dat <- gather(info.m, key = 'Year', value = 'Frequency', 2:26)
head(gat.dat)

i.1987 <- gat.dat %>% 
  filter(Year == 1987)

plot <- ggplot(i.1987, aes(x = index, y = Frequency)) + geom_bar(stat="identity") +
  theme_bw()

plot + theme(axis.line = element_line(),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             panel.border = element_blank(),
             panel.background = element_blank()
)


#### Barplot each year ###

par(mfrow = c(3,3))
for (i in years){
  freq <- gat.dat %>% 
    filter(Year == i)
  plot <- ggplot(freq, aes(x = index, y = Frequency)) + geom_bar(stat="identity") + 
    scale_y_continuous(expand = c(0,0)) +
    theme_bw()
  
  
  pdf(paste("D:/Stony Brook/Peconic Estuary/Peconic_Project/Barplots/Year_", i, ".pdf"), paper = "a4r")
  print(plot + theme(axis.line = element_line(),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     panel.border = element_blank(),
                     panel.background = element_blank()
  )) 
  dev.off()
}

#### Facet for each year ####
pgs <- ceiling(length(levels(factor(gat.dat$Year)))/9)

base <- ggplot(gat.dat, aes(x = index, y = Frequency)) + 
  geom_bar(stat="identity", color = "black") + 
  scale_y_continuous(breaks = round(seq(min(gat.dat$Frequency), max(gat.dat$Frequency), by = 2),1)) +
  scale_x_continuous(breaks = round(seq(min(gat.dat$index), max(gat.dat$index), by = 4),1)) +
  xlab("Station") + ylab("Frequency") + 
  theme_bw()

#pdf(paste("D:/Stony Brook/Peconic Estuary/Peconic_Project/Plots/Barplots/Yearly/All_Years.pdf"), paper = "a4r",
#    width = 10, height = 10)
for(i in seq_len(pgs)){
  print(base + facet_wrap_paginate(~Year, ncol = 3, nrow = 3, page = i) + 
          theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_blank(),
                panel.background = element_blank(),
                axis.text.x = element_text(vjust = 1, hjust=0.5, size = 6)
          ))
}
#dev.off()

#### Barplot each site  ####

index <- rownames(info.s)
sites.y <- cbind(index, info.s)
head(sites.y)

gat.y <- gather(sites.y, key = 'Quadrat', value = 'Frequency', 2:78)
head(gat.y)
gat.y[is.na(gat.y)] <- 0
gat.y$quad_f <- factor(gat.y$Quadrat, levels = 1:77)

par(mfrow = c(3,3))
for (i in stations){
  freq.y <- gat.y %>% 
    filter(Quadrat == i)
  plot.y <- ggplot(freq.y, aes(x = index, y = Frequency)) + geom_bar(stat="identity") + 
    scale_y_continuous(expand = c(0,0)) + xlab("Year") + ylab("Frequency") +
    theme_bw()
  
  #pdf(paste("D:/Stony Brook/Peconic Estuary/Peconic_Project/Plots/Barplots/Stations/Station_", i, ".pdf"), paper = "a4r")
  print(plot.y + theme(axis.line = element_line(),
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(),
                       panel.border = element_blank(),
                       panel.background = element_blank(),
                       axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)
  )) 
  #dev.off()
}

# Facet for each site #

face.y <- ggplot(gat.y, aes(x = index, y = Frequency)) + 
  geom_bar(stat="identity", color = "black") + 
  scale_y_continuous(breaks = round(seq(min(gat.y$Frequency), max(gat.y$Frequency), by = 2),1)) + 
  xlab("Year") + ylab("Frequency") +
  theme_bw()

pgs <- ceiling(length(levels(gat.y$quad_f))/9)

#pdf(paste("D:/Stony Brook/Peconic Estuary/Peconic_Project/Plots/Barplots/Stations/All_Stations.pdf"), paper = "a4r",
#    width = 10, height = 10)
for(i in seq_len(pgs)){
  print(face.y + facet_wrap_paginate(~quad_f, ncol = 3, nrow = 3, page = i) + 
          theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_blank(),
                panel.background = element_blank(),
                axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1, size = 8)
          ))
}
#dev.off()