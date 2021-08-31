
library(vegan)
library(tidyverse)
#library(mvabund)
#library(modeest)
library(ggplot2)
library(reshape2)
library(grid)
library(gridExtra)
#library(ggvegan)
library(ggrepel)
library(factoextra)
library(ggforce)

#### Total number of species caught per year per site ####

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

#### Generate data frame for trajectory Analaysis ####
merge.df.2 <- merge.df %>% arrange(Year) # arrange dataframe by year
merge.df.2 <- merge.df.2[,-1:-10]
merge.df.2 <- subset(merge.df.2, select = -c(Week, Month, Date, Season, StartLat, EndLat,
                           StartLong, EndLong, Grand.Total))
head(merge.df.2)

spp.names <- colnames(master.sheet[,2:42])
names(merge.df.2)[3:43] <- spp.names
merge.df.2[is.na(merge.df.2)] <- 0


#### Start aggregation loop for species only ####

x.sums <- NULL
xxx <- NULL
all.y.s <- NULL

for (i in years) {
  xxx <- merge.df.2 %>% 
    filter(Year == i)
  xxx <- xxx %>% arrange(Station)
  
  for (j in stations){
    xyx <- xxx %>% 
      filter(Station == j)
    xyx <- xyx[,-1:-2]
    col.sums <- colSums(xyx)
    x.sums <- rbind(x.sums, col.sums)
    yr.st <- cbind(as.data.frame(i), as.data.frame(j))
    all.y.s <- rbind(all.y.s, yr.st)
  }
}

all.y.s$Name <- unite(all.y.s, col=Name,i,j, sep = "_S")
all.y.s <- all.y.s[,-3]
colnames(all.y.s) <- c("Year", "Station")
x.sums <- as.data.frame(x.sums)
x.sums <- cbind(all.y.s, x.sums)
rownames(x.sums) <- 1:nrow(x.sums)

write.csv(x.sums, "D:/Stony Brook/Peconic Estuary/Peconic_Project/Data/Modified_Data/Total_per_site.csv")

#### Example
x.sums <- NULL

xxx <- merge.df.2 %>% 
  filter(Year == 1988)

xxx <- xxx %>% arrange(Station)
head(xxx)
xyx <- xxx %>%
  filter(Station == 1)
head(xyx)
temp.sum <- colSums(xyx)

x.sums <- rbind(x.sums, temp.sum)

for (i in stations){
  xyx <- xxx %>%
    filter(Station == i)
  xyx <- xyx[,-1:-2]
  col.sums <- colSums(xyx)
  x.sums <- rbind(x.sums, col.sums)
}

rownames(x.sums) <- stations
x.sums <- as.data.frame(x.sums)
x.sums$Year <- "1987"
x.sums$Station <- stations

x.temp <- x.sums
View(x.temp)


#### End of Example ####

#### Run aggregation loop for env. vars. and species ####
# Env. variables will be averaged per site #
# Species will be totaled per site #

df.try <- subset(merge.df, select = -c(id, SurfaceDO, Week, Month, Date, Season, StartLat, EndLat, 
                                       StartLong, EndLong, Grand.Total))
names(df.try)[11:51] <- spp.names
colnames(df.try[11:ncol(df.try)])
ncol(df.try)
spp.only <- data.frame(df.try[,11:ncol(df.try)])
spp.only[is.na(spp.only)] <- 0
head(spp.only)
total.row <- data.frame(rowSums(spp.only)) # no 0s
which(total.row == 0)
total.row

#write.csv(df.try, "D:/Stony Brook/Peconic Estuary/Peconic_Project/Data/Modified_Data/df_try.csv")

x.sums2 <- NULL
xxx2 <- NULL
all.y.s2 <- NULL
env.vars.mean <- NULL

for (i in years) {
  xxx2 <- df.try %>% 
    filter(Year == i)
  xxx2 <- xxx2 %>% arrange(Station)
  print(paste("Year ", i))
  
  for (j in unique(xxx2$Station)){
    sort.s <- xxx2 %>% arrange(Station)
    xyx2 <- sort.s %>% 
      filter(Station == j)
    print(paste("Station ", j))
    
    vars <- xyx2[, 1:8]
    vars.mean <- colMeans(vars)
    env.vars.mean <- rbind(env.vars.mean, vars.mean)
    
    xyx2 <- xyx2[,-1:-10]
    xyx2[is.na(xyx2)] <- 0
    col.sums2 <- colSums(xyx2)
    
    x.sums2 <- rbind(x.sums2, col.sums2)

    yr.st2 <- cbind(as.data.frame(i), as.data.frame(j))
    all.y.s2 <- rbind(all.y.s2, yr.st2)

  }
}

#all.y.s2$Name <- unite(all.y.s, col=Name,i,j, sep = "_S")
colnames(all.y.s2) <- c("Year", "Station")
x.sums2 <- as.data.frame(x.sums2)
x.sums2 <- cbind(all.y.s2, x.sums2)
rownames(x.sums2) <- 1:nrow(x.sums2)

env.vars.mean <- as.data.frame(env.vars.mean)
env.vars.mean <- cbind(all.y.s2, env.vars.mean)
rownames(env.vars.mean) <- 1:nrow(env.vars.mean)


write.csv(x.sums2, "D:/Stony Brook/Peconic Estuary/Peconic_Project/Data/Modified_Data/Total_per_site.csv")
write.csv(env.vars.mean, "D:/Stony Brook/Peconic Estuary/Peconic_Project/Data/Modified_Data/Env_vars_mean_per_site.csv")

#### Impute the missing means in env.vars ####
library(mice)
impute.vars <- mice(env.vars.mean[,3:ncol(env.vars.mean)], m = 10, method = "pmm", maxit = 100, seed = 500)
impute.vars <- complete(impute.vars, 6)
rownames(impute.vars) <- 1:nrow(impute.vars)
View(impute.vars)

impute.vars <- cbind(all.y.s2, impute.vars)

## Combine environmental data with abundance ##

total.mean.df <- cbind(impute.vars, x.sums2)
total.mean.df <- total.mean.df[,-11:-12]
total.mean.df$Year <- as.factor(total.mean.df$Year) # convert year as categorical variable

## Example ordination - need to remove columns with 0 total ##

dat <- total.mean.df %>%
  filter(Year == 1987)
head(dat)
dat <- subset(dat, select = -c(Year, Station))
env.dat.1987 <- dat[,1:8]
dat <- dat[,-1:-8]
dat[is.na(dat)] <- 0
dat <- dat[rowSums(dat[])>0,]


rda.n <- rda(dat)
rda.n
plot(rda.n)

pca.n <- prcomp(dat)
pca.n

fviz_pca_ind(pca.n,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE)     # Avoid text overlapping)

dca.s <- decostand(dat, method = "chi.square")
dca.s <- dca.s[, colSums(is.na(dca.s)) == 0] # chi squared computation require removal of columns with NAN
dca.n <- decorana(dca.s)
plot(dca.n)

temp <- envfit(dca.n, env.dat.1987, perm=999)

nmds.n <- metaMDS(dat)
plot(nmds.n)

#### Run nMDS ordination with just 77 points for each year ####
for (i in years){
  dat <- total.mean.df %>%
    filter(Year == i)
  head(dat)
  dat <- subset(dat, select = -c(Year, Station)) # can exclude year
  head(dat)
  
  env.dat.i <- dat[,1:8]
  dat <- dat[,-1:-8]
  dat[is.na(dat)] <- 0
  
  # Ordination Analysis #
  nMDS.i <- metaMDS(dat, distance = "bray", autotransform = F)
  
  # Species ordination
  spp.fit <- envfit(nMDS.i, dat, permutations = 999) # fits species vectors
  # Site ordination
  site.scrs <- as.data.frame(scores(nMDS.i, display = "sites"))
  head(site.scrs)
  
  ## site.scrs <- cbind(site.scrs, Management = dune.env$Management) #add grouping variable "Management" to dataframe
  ## site.scrs <- cbind(site.scrs, Landuse = dune.env$Use) #add grouping variable of cluster grouping to dataframe
  ## site.scrs <- cbind(site.scrs, Site = rownames(site.scrs)) #add site names as variable if you want to display on plot
  
  spp.scrs <- as.data.frame(scores(spp.fit, display = "vectors")) #save species intrinsic values into dataframe
  spp.scrs <- cbind(spp.scrs, Species = rownames(spp.scrs)) # add species name to vectors
  spp.scrs <- cbind(spp.scrs, pval = spp.fit$vectors$pvals) #add pvalues to dataframe so you can select species which are significant
  sig.spp.scrs <- subset(spp.scrs, pval<=0.05)
  head(spp.scrs)
  
  #### Environmental Variables ####
  
  # Fit nMDS to environmental data #
  
  i.envfit <- envfit(nMDS.i, env.dat.i, permutations = 999) # this fits environmental vectors
  
  env.i.scores <- as.data.frame(scores(i.envfit, display = "vectors")) #extracts relevant scores from envifit
  env.i.scores <- cbind(env.i.scores, env.variables = rownames(env.i.scores)) #and then gives them their names
  env.i.scores <- cbind(env.i.scores, pval = i.envfit$vectors$pvals) # add pvalues to dataframe
  sig.env.scrs <- subset(env.i.scores, pval<=0.05) #subset data to show variables significant at 0.05
  head(sig.env.scrs)
  
  #### GGPLOT ####
  # Plot sites
  nmds.plot <- ggplot(site.scrs, aes(x=NMDS1, y=NMDS2))+ #sets up the plot
    geom_text(aes(NMDS1, NMDS2, size = 1), color = "grey", label= row.names(site.scrs), size = 2) + #adds site points to plot, shape determined by Landuse, colour determined by Management
    coord_fixed()+
    theme_classic()+ 
    theme(panel.background = element_rect(fill = NA, colour = "black", size = 1, linetype = "solid"))+
    #labs() + # add legend labels for Management and Landuse
    theme(legend.position = "right", legend.text = element_text(size = 12), 
          legend.title = element_text(size = 12), axis.text = element_text(size = 10),
          axis.line = element_line(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank()
    ) 
  
  
  pdf(paste("D:/Stony Brook/Peconic Estuary/PeconicEstuary/Ordination/nMDS_mean/Bi-Tri_plots/Chisq/Year_", i, ".pdf"), paper = "a4r",
      width = 10, height = 10, onefile = TRUE)
  
  print(nmds.plot +
          geom_segment(data = sig.spp.scrs, aes(x = 0, xend=NMDS1, y=0, yend=NMDS2), 
                       arrow = arrow(length = unit(0.25, "cm")), colour = "red", lwd=0.25) + #add vector arrows of significant species
          ggrepel::geom_text_repel(data = sig.spp.scrs, aes(x=NMDS1, y=NMDS2, label = Species), 
                                   color = "red", cex = 4, direction = "both", segment.size = 0.5) #add labels for species, use ggrepel::geom_text_repel so that labels do not overlap
  )
  
  print(nmds.plot +
          geom_segment(data = env.i.scores, aes(x = 0, xend=NMDS1, y=0, yend=NMDS2), 
                       arrow = arrow(length = unit(0.25, "cm")), colour = "black", lwd=0.25) + #add vector arrows of significant species
          ggrepel::geom_text_repel(data = env.i.scores, aes(x=NMDS1, y=NMDS2, label = env.variables), max.overlaps = 20,
                                   color = "black", cex = 3, direction = "both", segment.size = 0.5))
  dev.off()
}

#### Run DCA ordination ####  

for (i in years){
  dat <- total.mean.df %>%
    filter(Year == i)
  head(dat)
  dat <- subset(dat, select = -c(Year, Station)) # can exclude year
  head(dat)
  
  env.dat.i <- dat[,1:8]
  dat <- dat[,-1:-8]
  dat[is.na(dat)] <- 0

  dca.s <- decostand(dat, method = "chi.square")
  dca.s <- dca.s[, colSums(is.na(dca.s)) == 0] # chi squared computation require removal of columns with NAN
  dca.n <- decorana(dca.s)
  plot(dca.n, main= paste("Year", i))
  
  # Species ordination
  spp.fit <- envfit(nMDS.i, dat, permutations = 999) # fits species vectors
  # Site ordination
  site.scrs <- as.data.frame(scores(nMDS.i, display = "sites"))
  head(site.scrs)
  
  ## site.scrs <- cbind(site.scrs, Management = dune.env$Management) #add grouping variable "Management" to dataframe
  ## site.scrs <- cbind(site.scrs, Landuse = dune.env$Use) #add grouping variable of cluster grouping to dataframe
  ## site.scrs <- cbind(site.scrs, Site = rownames(site.scrs)) #add site names as variable if you want to display on plot
  
  spp.scrs <- as.data.frame(scores(spp.fit, display = "vectors")) #save species intrinsic values into dataframe
  spp.scrs <- cbind(spp.scrs, Species = rownames(spp.scrs)) # add species name to vectors
  spp.scrs <- cbind(spp.scrs, pval = spp.fit$vectors$pvals) #add pvalues to dataframe so you can select species which are significant
  sig.spp.scrs <- subset(spp.scrs, pval<=0.05)
  head(spp.scrs)
  
  #### Environmental Variables ####
  
  # Fit nMDS to environmental data #
  
  i.envfit <- envfit(nMDS.i, env.dat.i, permutations = 999) # this fits environmental vectors
  
  env.i.scores <- as.data.frame(scores(i.envfit, display = "vectors")) #extracts relevant scores from envifit
  env.i.scores <- cbind(env.i.scores, env.variables = rownames(env.i.scores)) #and then gives them their names
  env.i.scores <- cbind(env.i.scores, pval = i.envfit$vectors$pvals) # add pvalues to dataframe
  sig.env.scrs <- subset(env.i.scores, pval<=0.05) #subset data to show variables significant at 0.05
  head(sig.env.scrs)
  
  #### GGPLOT ####
  # Plot sites
  nmds.plot <- ggplot(site.scrs, aes(x=NMDS1, y=NMDS2))+ #sets up the plot
    geom_text(aes(NMDS1, NMDS2, size = 1), color = "grey", label= row.names(site.scrs), size = 2) + #adds site points to plot, shape determined by Landuse, colour determined by Management
    coord_fixed()+
    theme_classic()+ 
    theme(panel.background = element_rect(fill = NA, colour = "black", size = 1, linetype = "solid"))+
    #labs() + # add legend labels for Management and Landuse
    theme(legend.position = "right", legend.text = element_text(size = 12), 
          legend.title = element_text(size = 12), axis.text = element_text(size = 10),
          axis.line = element_line(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank()
    ) 
  
  
  pdf(paste("D:/Stony Brook/Peconic Estuary/PeconicEstuary/Ordination/nMDS_mean/Bi-Tri_plots/Chisq/Year_", i, ".pdf"), paper = "a4r",
      width = 10, height = 10, onefile = TRUE)
  
  print(nmds.plot +
          geom_segment(data = sig.spp.scrs, aes(x = 0, xend=NMDS1, y=0, yend=NMDS2), 
                       arrow = arrow(length = unit(0.25, "cm")), colour = "red", lwd=0.25) + #add vector arrows of significant species
          ggrepel::geom_text_repel(data = sig.spp.scrs, aes(x=NMDS1, y=NMDS2, label = Species), 
                                   color = "red", cex = 4, direction = "both", segment.size = 0.5) #add labels for species, use ggrepel::geom_text_repel so that labels do not overlap
  )
  
  print(nmds.plot +
          geom_segment(data = env.i.scores, aes(x = 0, xend=NMDS1, y=0, yend=NMDS2), 
                       arrow = arrow(length = unit(0.25, "cm")), colour = "black", lwd=0.25) + #add vector arrows of significant species
          ggrepel::geom_text_repel(data = env.i.scores, aes(x=NMDS1, y=NMDS2, label = env.variables), max.overlaps = 20,
                                   color = "black", cex = 3, direction = "both", segment.size = 0.5))
  dev.off()
}

}
