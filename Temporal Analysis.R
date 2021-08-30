#### Temporal Analysis ####

library(devtools)
install_github("james-thorson/VAST", INSTALL_opts="--no-staged-install")
install.packages("codyn")

library(VAST)
library(codyn)

#### Using *codyn* ####

## Species turnover
## Total turnover

data("collins08")
View(collins08)

# `turnover` requires a data frame with columns for species, time and abundance, 
# and includes an optional argument to specify a column for spatial replicates.

KNZ_turnover <- turnover(df = collins08, 
                         time.var = "year", 
                         species.var = "species", 
                         abundance.var = "abundance", 
                         replicate.var = "replicate")
View(KNZ_turnover)

# Specifying metric="appearance" will return the proportion of species that appeared 
# relative to the total number of species observed in both time points.

KNZ_appearance <- turnover(df = collins08, 
                           time.var = "year",
                           species.var = "species",
                           abundance.var = "abundance",
                           replicate.var = "replicate",
                           metric = "appearance")
View(KNZ_appearance)

# Similarly, specifying metric="disappearance" will return the proportion of species that disappeared 
# relative to the total number of species observed in both time points.

KNZ_disappearance <- turnover(df = collins08,
                              time.var = "year",
                              species.var = "species",
                              abundance.var = "abundance",
                              replicate.var = "replicate",
                              metric = "disappearance")
View(KNZ_disappearance)

plot <- ggplot(KNZ_turnover, aes(x=year, y=total)) + 
  geom_line() +
  geom_point() +
  theme_bw() + theme(legend.key=element_blank())

plot = ggplot(data = KNZ_turnover, aes(x = year, y = total, col = replicate)) + geom_line(lwd = 2) +
  theme(legend.key=element_blank()) + 
  xlab('Year') +
  ylab('Total Abundance') + 
  geom_point ()

plot

plot + theme(axis.line = element_line(),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             panel.border = element_blank(),
             panel.background = element_blank()
)
