plotTrajDiamDist<-function(cli = 7) {
  l = colnames(avoca_strat[[1]])
  ncl = 14
  m197072= avoca_strat[avoca_surveys==1][[cli]]["NOTCLI",2:ncl]
  m197072[m197072<1] = NA
  m1974 = avoca_strat[avoca_surveys==2][[cli]]["NOTCLI",2:ncl]
  m1974[m1974<1] = NA
  m1978 = avoca_strat[avoca_surveys==3][[cli]]["NOTCLI",2:ncl]
  m1978[m1978<1] = NA
  m1983 = avoca_strat[avoca_surveys==4][[cli]]["NOTCLI",2:ncl]
  m1983[m1983<1] = NA
  m1987 = avoca_strat[avoca_surveys==5][[cli]]["NOTCLI",2:ncl]
  m1987[m1987<1] = NA
  m1993 = avoca_strat[avoca_surveys==6][[cli]]["NOTCLI",2:ncl]
  m1993[m1993<1] = NA
  m1999 = avoca_strat[avoca_surveys==7][[cli]]["NOTCLI",2:ncl]
  m1999[m1999<1] = NA
  m2004 = avoca_strat[avoca_surveys==8][[cli]]["NOTCLI",2:ncl]
  m2004[m2004<1] = NA
  m2009 = avoca_strat[avoca_surveys==9][[cli]]["NOTCLI",2:ncl]
  m2009[m2009<1] = NA
  
  
  plot(m197072, type="l", ylim=c(1,200), log="y",
       xlab="", ylab="Number of individuals (log)", main=paste0("Trajectory ",cli), 
       axes=FALSE, col=gray(0.8), lwd=2)
  axis(2, las=2)
  axis(1, at=1:(ncl-1), labels=l[2:ncl], las=2)
  lines(m1974, col=gray(0.7), lwd=2)
  lines(m1978, col=gray(0.6), lwd=2)
  lines(m1983, col=gray(0.5), lwd=2)
  lines(m1987, col=gray(0.4), lwd=2)
  lines(m1993, col=gray(0.3), lwd=2)
  lines(m1999, col=gray(0.2), lwd=2)
  lines(m2004, col=gray(0.1), lwd=2)
  lines(m2009, col=gray(0), lwd=2)
  legend("topright", bty="n", lwd=2,col=gray(seq(0.8,0, by=-0.1)), legend=c("1970/72","1974","1978","1983", "1987", "1993","1999","2004","2009"))
}
