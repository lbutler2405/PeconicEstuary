#### Data Preprocessing ####
install.packages("modeest")
#install.packages("ggvegan")
install.packages("ggrepel")
install.packages("factoextra")
install.packages("ggforce")

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

## Station 1
st.1 <- merge.df %>% 
  filter(Station == 1) 

head(st.1)
nrow(st.1)

View(st.1)
unique(st.1$Year)
table.1 <- table(st.1$Year)
table.1

## Station 2

st.2 <- merge.df %>% 
  filter(Station == 2) 

head(st.2)
nrow(st.2)

View(st.2)
unique(st.2$Year)
table.2 <- table(st.2$Year)
table.2


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

## Find median, min, max and mean for each row

min.r <- data.frame(apply(info.df, 1, min))
max.r <- data.frame(apply(info.df, 1, max))
mean.r <- data.frame(apply(info.df, 1, mean))
median.r <- data.frame(apply(info.df, 1, median))

calc.rr <- cbind(min.r, max.r, mean.r, median.r)
colnames(calc.rr) <- c("Min", "Max", "Mean", "Median")
calc.rr[79:82, ] <- NA


info <- cbind(info, calc.rr)
write.csv(info, "info_spp_yr.csv")

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

## plot frequency for each site ##

#melt data frame into long format
index <- 1:nrow(info.df)
info.m <- cbind(index, info.df)
#data.m <- melt(info.m, id.vars = 'index', variable.name = 'series')

# plot all years - too messy

ggplot(data.melt, aes(index, value)) +
  geom_line(aes(colour = series)) + 
  theme_bw()

# Plot CCA for each year

yr.cca <- cca(info.df)
par(mfrow = c(1,2))
plot(yr.cca, display = "sites") # site 59 produces issues due to large number of 0s (i.e. they did not go to site 59)
plot(yr.cca, display = "species")

# removing site 59

yr.cca.2 <- cca(info.df[-59,])
par(mfrow = c(1,2))
plot(yr.cca.2, display = "sites", main = "Quadrats")
plot(yr.cca.2, display = "species", main = "Years")

#### Number of times sample was taken per season per year ####

years

info.y <- NULL

for (i in unique(years)){
  st.i <- merge.df %>% 
    filter(Year == i)
  table. <- table(st.i$Month)
  info.y <- dplyr::bind_rows(info.y, table.)
}

head(info.y)
rownames(info.y) <- years

#### Number of times sample was taken per week per site ####

info.w <- NULL

for (i in unique(week)){
  st.i <- merge.df %>% 
    filter(Week == i)
  table. <- table(st.i$Station)
  info.w <- dplyr::bind_rows(info.w, table.)
}

head(info.w)
rownames(info.w) <- week

colnames.w <- colnames(info.w)
colnames.w
temp3 <- t(info.w)
temp3 <- as.data.frame(temp3, col.names = col.names(temp3))
head(temp3)
temp3$index <- as.numeric(row.names(temp3))
head(temp3)

temp3 <- temp3[order(temp3$index), ]
head(temp3)
temp3 <- temp3[,-29]
temp3 <- t(temp3)
head(temp3)
info.w <- temp3

info.w[is.na(info.w)] <- 0

# Calculate mode for each row
mode.w <- NULL

for (i in 1:nrow(info.w)){
  mode.i <- getmode(info.w[i,])
  mode.w <- rbind(mode.w, mode.i)
}

mode.w
info.w <- cbind(info.w, mode.w)
info.w <- as.data.frame(info.w)
names(info.w)[78] <- 'Mode'


info.w <- t(info.w)
# remove mode
info.w <- info.w[-78, ]

## CCA on weeks vs. quadrats
week.cca <- cca(info.w)

par(mfrow=c(1,2))
plot(week.cca, display = "species") # weeks
plot(week.cca, display = "sites") # quadrats

# remove site 59
info.w.2 <- info.w[-59,] 
week.cca.2 <- cca(info.w.2)

plot(week.cca.2, display = "species", main = "Weeks") # weeks
plot(week.cca.2, display = "sites", main = "Quadrats") # quadrats


#### save data sorted by week for each year in a separate .csv file ####

for (i in years){
  st.i <- merge.df %>% 
    filter(Year == i)
  st.i <- st.i %>% arrange(Week)
  write.csv(st.i, paste('Example_Data/Data_of_Year_', i, '.csv'))
}

#####
#### save data sorted by station for each year in a separate .csv file ####

for (i in years){
  st.s <- merge.df %>% 
    filter(Year == i)
  st.s <- st.s %>% arrange(Station)
  write.csv(st.s, paste('Example_Data/sorted_station/Data_of_Year_', i, '.csv'))
}

#### Generate data frame for trajectory Analaysis ####
merge.df.2 <- merge.df %>% arrange(Year) # arrange dataframe by year
merge.df.2 <- merge.df.2[,-17:-20]

merge.df.2 <- merge.df.2[, -1:-10]

dates <- merge.df.2[,1:6]
head(dates)
merge.df.2 <- merge.df.2[, -1:-6]

merge.df.2 <- merge.df.2[,-42] # drop grand total
merge.df.3 <- merge.df.2
merge.df.3[is.na(merge.df.3)] <- 0

# run cca on species data #

spp.cca <- cca(merge.df.3)
par(mfrow = c(1,2))
plot(spp.cca, display = "species")
plot(spp.cca, display = "sites")

# run DCA - CCA provided wedges #
spp.stan <- decostand(merge.df.3, method = 'hellinger') # standardize
spp.dca <- decorana(spp.stan)
head(spp.dca)
head(spp.dca$rproj)
head(spp.dca$cproj)

plot(spp.dca, display = "sites")
plot(spp.dca, display = "species")

# nMDS #

spp.nmds <- metaMDS(merge.df.3)

# RDA #
spp.rda <- rda(merge.df.3)
par(mfrow=c(1,2))
plot(spp.rda, display = "species")
plot(spp.rda, display = "sites")

# Using hellinger standardized data
spp.rda.2 <- rda(spp.stan)
par(mfrow=c(1,2))
plot(spp.rda.2, display = "species")
plot(spp.rda.2, display = "sites")

#### PCA #### Dataset is too large

spp.pca <- prcomp(merge.df.3)
fviz_pca_ind(spp.pca,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE)     # Avoid text overlapping)

## To get values from PCA ##
# Eigenvalues

eig.val <- get_eigenvalue(res.pca)
eig.val

# Results for Variables
res.var <- get_pca_var(res.pca)
res.var$coord          # Coordinates
res.var$contrib        # Contributions to the PCs
res.var$cos2           # Quality of representation 
# Results for individuals
res.ind <- get_pca_ind(res.pca)
res.ind$coord          # Coordinates
res.ind$contrib        # Contributions to the PCs
res.ind$cos2           # Quality of representation 

#### PCA for each year ####

spp.dat.yr <- cbind(dates$Year, merge.df.3)
names(spp.dat.yr)[1] <- "Year"
spp.names <- colnames(master.sheet[,2:42])
names(spp.dat.yr)[2:42] <- spp.names
head(spp.dat.yr)

par(mfrow=c(1,2))

for (i in years){
  st.i <- spp.dat.yr %>% 
    filter(Year == i)
  st.i <- st.i[,-1]
  pca <- prcomp(st.i)
  
  # plot individuals
  fviz_pca_ind(pca,
               col.ind = "cos2", # Color by the quality of representation
               gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
               repel = TRUE)     # Avoid text overlapping
  
  # plot variables (species)
  
  fviz_pca_var(pca,
               col.var = "contrib", # Color by contributions to the PC
               gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
               repel = TRUE)     # Avoid text overlapping
  
  eig.val <- get_eigenvalue(pca)
  res.var <- get_pca_var(pca) # results for variables
  coords_var <- res.var$coord  
  
  res.ind <- get_pca_ind(pca) # results for individuals
  coords_ind <- res.ind$coord          # Coordinates
  
  write.csv(eig.val, paste('Example_Data/CCA_results/PCA_Eigen_', i, '.csv'))
  write.csv(coords_var, paste('Example_Data/CCA_results/PCA_variables_', i, '.csv'))
  write.csv(coords_ind, paste('Example_Data/CCA_results/PCA_individuals_', i, '.csv'))
}

#### CCA for each year ####
par(mfrow = c(1,2))
for (i in years){
  st.i <- spp.dat.yr %>% 
    filter(Year == i)
  st.i <- st.i[,-1]
  cca <- cca(st.i)

  plot(cca, display = "sites", main =  paste("Site for ",i))
  plot(cca, display = "species", main = paste("Species for ", i))     
}


#### DCA for each year #### - MAYBE
par(mfrow = c(2,2), mar = c(2, 2, 1.2, 1.2))
for (i in years){
  st.i <- spp.dat.yr %>% 
    filter(Year == i)
  st.i <- st.i[,-1]
  decos <- decostand(st.i, method = "hellinger")
  dca <- decorana(decos)
  par(mfrow=c(2,2), mar = c(2, 2, 1.2, 1.2))

  pdf(paste("D:/Stony Brook/Peconic Estuary/Peconic_Project/DCA_old/Year_", i, ".pdf"), paper = "a4r")
  par(mfrow=c(2,2), mar = c(2, 2, 1.2, 1.2))
  print(plot(dca, type = "n", display = "sites", main =  paste("Site for ",i)))
  print(text(dca, display = "sites", cex = 0.5, pos=2))
  print(plot(dca, display = "species", cex = 0.5, main = paste("Species for ", i)))
  print(plot(dca, type = "n", main = paste("Biplot for", i)))
  print(text(dca, display = "sites", cex = 0.5, pos = 2))
  print(text(dca, display = "species", cex = 0.5, pos = 1, col = "red"))  
  dev.off()
}

#### RDA for each year #### - NO

par(mfrow = c(1,2))
for (i in years){
  st.i <- spp.dat.yr %>% 
    filter(Year == i)
  st.i <- st.i[,-1]
  rda <- rda(st.i)
  plot(rda, display = "sites", main =  paste("Site for ",i))
  plot(rda, display = "species", main = paste("Species for ", i))     
}

#### NMDS for each year ####

par(mfrow = c(1,2))
for (i in years){
  st.i <- spp.dat.yr %>% 
    filter(Year == i)
  st.i <- st.i[,-1]
  nMDS <- metaMDS(st.i)
  
  pdf(paste("D:/Stony Brook/Peconic Estuary/Peconic_Project/nMDS_old/Year_", i, ".pdf"), paper = "a4r")
  par(mfrow=c(1,2), mar = c(0.1, 0.2, 0.2, 0.1))
  print(plot(nMDS, display = "sites", main =  paste("Site for ",i), type = "t"))
  print(plot(nMDS, display = "species", main = paste("Species for ", i), type = "t")) 
  dev.off()
}

#### NMDS Site Scores ####
# Example

i.1987 <- spp.dat.yr %>%
  filter(Year == 1987)
head(i.1987)
i.1987 <- i.1987[,-1]

nMDS.1987 <- metaMDS(i.1987)
spp.fit <- envfit(nMDS.1987, i.1987, permutations = 999) # fits species vectors

site.scrs <- as.data.frame(scores(nMDS.1987, display = "sites"))
head(site.scrs)

## site.scrs <- cbind(site.scrs, Management = dune.env$Management) #add grouping variable "Management" to dataframe
## site.scrs <- cbind(site.scrs, Landuse = dune.env$Use) #add grouping variable of cluster grouping to dataframe
## site.scrs <- cbind(site.scrs, Site = rownames(site.scrs)) #add site names as variable if you want to display on plot

spp.scrs <- as.data.frame(scores(spp.fit, display = "vectors")) #save species intrinsic values into dataframe
spp.scrs <- cbind(spp.scrs, Species = rownames(spp.scrs)) # add species name to vectors
spp.scrs <- cbind(spp.scrs, pval = spp.fit$vectors$pvals) #add pvalues to dataframe so you can select species which are significant
head(spp.scrs)
sig.spp.scrs <- subset(spp.scrs, pval<=0.05) #subset data to show species significant at 0.05
head(sig.spp.scrs)

#### NMDS Species Scores ####



#### Add centroids - centroids are years ####
x <- results.dca$DCA1
y <- results.dca$DCA2
class <- quads.temp$Comm
df <- data.frame(class, x, y)

# create a dataframe of sequential years and sequential week numbers


main.week <- NULL

for (i in years) {
  yr <- merge.df.2 %>% 
    filter (Year == i)
  yr <- yr %>% arrange(Week)
  main.week <- rbind(main.week, yr)
}

main.st <- NULL

for (i in years) {
  wk <- merge.df.2 %>% 
    filter (Year == i)
  wk <- wk %>% arrange(Station)
  main.st <- rbind(main.st, wk)
}

##### Barplots - frequency of each site for each year ####

### Sample ####

index <- 1:nrow(info.df)
info.m <- cbind(index, info.df)

gat.dat <- gather(info.m, key = 'Year', value = 'Frequency', 2:26)
head(gat.dat)

i.f <- gat.dat %>% 
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

pdf(paste("D:/Stony Brook/Peconic Estuary/Peconic_Project/Plots/Barplots/Yearly/All_Years.pdf"), paper = "a4r",
    width = 10, height = 10)
for(i in seq_len(pgs)){
print(base + facet_wrap_paginate(~Year, ncol = 3, nrow = 3, page = i) + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(vjust = 1, hjust=0.5, size = 6)
        ))
}
dev.off()

#### Barplot each site  ####

index <- rownames(info.s)
sites.y <- cbind(index, info.s)

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
  
  pdf(paste("D:/Stony Brook/Peconic Estuary/Peconic_Project/Plots/Barplots/Stations/Station_", i, ".pdf"), paper = "a4r")
  print(plot.y + theme(axis.line = element_line(),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     panel.border = element_blank(),
                     panel.background = element_blank(),
                     axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)
  )) 
  dev.off()
}

# Facet for each site #

face.y <- ggplot(gat.y, aes(x = index, y = Frequency)) + 
  geom_bar(stat="identity", color = "black") + 
  scale_y_continuous(breaks = round(seq(min(gat.y$Frequency), max(gat.y$Frequency), by = 2),1)) + 
  xlab("Year") + ylab("Frequency") +
  theme_bw()

pgs <- ceiling(length(levels(gat.y$quad_f))/9)

pdf(paste("D:/Stony Brook/Peconic Estuary/Peconic_Project/Plots/Barplots/Stations/All_Stations.pdf"), paper = "a4r",
    width = 10, height = 10)
for(i in seq_len(pgs)){
  print(face.y + facet_wrap_paginate(~quad_f, ncol = 3, nrow = 3, page = i) + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1, size = 8)
        ))
}
dev.off()


#### Environmental Data ####

colSums(is.na(env.data)) # Number of missing values per column

env.data <- env.data[,-17] # Removes the 'Comments' columns

env.data <- subset(env.data, select = -c(SURF_DO100,BOT_DO100)) # removes the columns for DO100 surface and bottom