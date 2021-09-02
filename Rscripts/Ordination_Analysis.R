#### Ordination Analysis and Plotting ####

#### Data Processing ####

#### 1. Example ####

# Data Processing #

# .csv file containing % abundance
master.sheet <- read.csv("species_master_spreadsheet.csv")
head(master.sheet)

# .csv file containing real abundances
raw.spp <- read.csv("all_spp_raw.csv")
head(raw.spp)

merge.df <- merge(master.sheet, raw.spp, by = "id")
head(merge.df)

merge.df <- merge.df[,-2:-42]
merge.df.2 <- merge.df %>% arrange(Year) # arrange dataframe by year
#merge.df.2 <- merge.df.2[,-17:-20]


# 1. Example

i.1987 <- merge.df.2 %>%
  filter(Year == 1987)
head(i.1987)
i.1987 <- i.1987[,-1]
i.1987 <- subset(i.1987, select = -c(Year, Week, Month, Date, Station, Season, StartLat, EndLat,
                                       StartLong, EndLong, Grand.Total))
head(i.1987)

spp.names <- colnames(master.sheet[,2:42])
names(i.1987)[10:50] <- spp.names

env.dat.1987 <- i.1987[,1:9]

i.1987 <- i.1987[,-1:-9]
i.1987[is.na(i.1987)] <- 0

# Ordination Analysis #

nMDS.1987 <- metaMDS(i.1987, distance = "bray", autotransform = F)


#### Plot nMDS ####
par(mfrow = c(2,2))
plot(nMDS.1987, display = "sites", main =  paste("Site for 1987"), type = "t")
plot(nMDS.1987, display = "species", main = paste("Species for 1987"), type = "t") 
plot(nMDS.1987, main = paste("BiPlot for 1987"), type = "t") 


#### Extract Values from nMDS ####
spp.fit.1987 <- envfit(nMDS.1987, i.1987, permutations = 999) # fits species vectors

site.scrs.1987 <- as.data.frame(scores(nMDS.1987, display = "sites"))
head(site.scrs.1987)

## site.scrs <- cbind(site.scrs, Management = dune.env$Management) #add grouping variable "Management" to dataframe
## site.scrs <- cbind(site.scrs, Landuse = dune.env$Use) #add grouping variable of cluster grouping to dataframe
## site.scrs <- cbind(site.scrs, Site = rownames(site.scrs)) #add site names as variable if you want to display on plot

spp.scrs.1987 <- as.data.frame(scores(spp.fit.1987, display = "vectors")) #save species intrinsic values into dataframe
spp.scrs.1987 <- cbind(spp.scrs.1987, Species = rownames(spp.scrs.1987)) # add species name to vectors
spp.scrs.1987 <- cbind(spp.scrs.1987, pval = spp.fit.1987$vectors$pvals) #add pvalues to dataframe so you can select species which are significant
head(spp.scrs.1987)
sig.spp.scrs.1987 <- subset(spp.scrs.1987, pval<=0.05) #subset data to show species significant at 0.05
head(sig.spp.scrs.1987)


#### Environmental Variables ####

# Fit nMDS using 1987 as example #

i1987.envfit <- envfit(nMDS.1987, env.dat.1987, permutations = 999) # this fits environmental vectors

env.1987.scores <- as.data.frame(scores(i1987.envfit, display = "vectors")) #extracts relevant scores from envifit
env.1987.scores <- cbind(env.1987.scores, env.variables = rownames(env.1987.scores)) #and then gives them their names
env.1987.scores <- cbind(env.1987.scores, pval = i1987.envfit$vectors$pvals) # add pvalues to dataframe

sig.env.scrs.1987 <- subset(env.1987.scores, pval<=0.05) #subset data to show variables significant at 0.05
head(env.1987.scores)

#### GGPLOT ####
# Plot sites
nmds.plot <- ggplot(site.scrs.1987, aes(x=NMDS1, y=NMDS2))+ #sets up the plot
  geom_text(aes(NMDS1, NMDS2, size = 1), color = "darkgrey", label= row.names(site.scrs.1987), size = 3) + #adds site points to plot, shape determined by Landuse, colour determined by Management
               coord_fixed()+
               theme_bw()+ 
               theme(panel.background = element_rect(fill = NA, colour = "black", size = 1, linetype = "solid"))+
               #labs() + # add legend labels for Management and Landuse
               theme(legend.position = "right", legend.text = element_text(size = 12), 
                     legend.title = element_text(size = 12), axis.text = element_text(size = 10),
                     axis.line = element_line(),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     panel.border = element_blank(),
                     panel.background = element_blank()
                     ) # add legend at right of plot

nmds.plot

# Add species

nmds.plot +
  geom_segment(data = spp.scrs.1987, aes(x = 0, xend=NMDS1, y=0, yend=NMDS2), 
               arrow = arrow(length = unit(0.25, "cm")), colour = "red", lwd=0.25) + #add vector arrows of significant species
  ggrepel::geom_text_repel(data = spp.scrs.1987, aes(x=NMDS1, y=NMDS2, label = Species), max.overlaps = 20,
                           color = "red", cex = 3, direction = "both", segment.size = 0.5) #add labels for species, use ggrepel::geom_text_repel so that labels do not overlap

#Add env.vars #
nmds.plot +
  geom_segment(data = env.1987.scores, aes(x = 0, xend=NMDS1, y=0, yend=NMDS2), 
               arrow = arrow(length = unit(0.25, "cm")), 
               colour = "black", lwd=0.3) + #add vector arrows of env variables
  ggrepel::geom_text_repel(data = env.1987.scores, aes(x=NMDS1, y=NMDS2, label = env.variables), cex = 3, 
                           direction = "both", segment.size = 0.25)#add labels for env variables
#########
# Add significant species only
nmds.plot +
  geom_segment(data = sig.spp.scrs.1987, aes(x = 0, xend=NMDS1, y=0, yend=NMDS2), 
               arrow = arrow(length = unit(0.25, "cm")), colour = "red", lwd=0.25) + #add vector arrows of significant species
  ggrepel::geom_text_repel(data = sig.spp.scrs.1987, aes(x=NMDS1, y=NMDS2, label = Species), max.overlaps = 20,
                           color = "red", cex = 3, direction = "both", segment.size = 0.5) #add labels for species, use ggrepel::geom_text_repel so that labels do not overlap


# Add significant Env. Vars. only

nmds.plot +
  geom_segment(data = sig.env.scrs.1987, aes(x = 0, xend=NMDS1, y=0, yend=NMDS2), 
               arrow = arrow(length = unit(0.25, "cm")), 
               colour = "black", lwd=0.3) + #add vector arrows of env variables
  ggrepel::geom_text_repel(data = sig.env.scrs.1987, aes(x=NMDS1, y=NMDS2, label = env.variables), cex = 3, 
                           direction = "both", segment.size = 0.25)#add labels for env variables

#### 2. Loop for each Year ####
for (i in years){
  dat <- x.sums %>%
    filter(Year == i)
  head(dat)
  dat <- dat[,-1]
  dat <- subset(dat, select = -c(Year, Week, Month, Date, Station, Season, StartLat, EndLat,
                                     StartLong, EndLong, Grand.Total))
  head(dat)
  
  spp.names <- colnames(master.sheet[,2:42])
  names(dat)[10:50] <- spp.names
  
  env.dat.i <- dat[,1:9]
  dat <- dat[,-1:-9]
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
  
  
  pdf(paste("D:/Stony Brook/Peconic Estuary/PeconicEstuary/Ordination/nMDS_old/Bi-Tri_plots/Year_", i, ".pdf"), paper = "a4r",
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



##### Simple Ordination plots ####

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




