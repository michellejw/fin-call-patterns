# Create plots for Call patterns paper
#library(ggplot2)
library(scales)
library(viridis)


# Load call patterns CSV
setwd('~/Dropbox/Research/CODE/CountPatterns/fin-call-patterns/SEQ_CODE/EXPLORE_DATA/CALLPATTERNS_STATS/CALLPATTERNS_STATS/')
dat = read.csv("../../ALL_seq_summary.csv", header = TRUE)
locs = read.csv("../../stationlist.csv", header = TRUE)

# Set up station colors
cpal_sta = c('#8dd3c7','#ffffb3','#bebada','#fb8072','#80b1d3')
cpal_sta = viridis(5)


# Convert datevec to R date format
dat$datevec <- as.Date(dat$datevec, "%Y-%m-%d")
dat$month <- as.numeric(format(dat$datevec,'%m'))
dat$year <- as.numeric(format(dat$datevec,'%Y'))

#mod_f5 <- lm(year ~ freq, data=dat)

#---------- Figure 5: 2 panels, Decadal F/IPI -----------

# Subset dataframe to include only Decadal analysis stations
decadedf = subset(dat,(station=='AX' | station=='KENE' | station=='KEMF') 
                  & peakcounts > 50)
# Further subset to include only notes in the low IPI/high freq group
decadedf_sub = subset(decadedf,(freq<21)&(ipi>22)) 
# Regression through decadal frequency
mod_f5_freq <- lm(freq~year, data=decadedf_sub)
mod_f5_ipi <- lm(ipi~year, data=decadedf_sub)


flims = c(17,25)
ipilims = c(10,40)

pdf(file="test.pdf",width=6,height=4)
fig5 = layout(matrix(c(1,2,1,3),2,2,byrow = TRUE))
par(oma=c(3,3,3,3), mar=c(0,.5,0,.5),mgp=c(1.4,.4,0))

# Panel 1
plot(decadedf$freq,decadedf$ipi,
     cex=decadedf$propcounts*.3,
     col=cpal_sta[decadedf$netvec],axes=FALSE,
     xlim=flims,ylim=ipilims)
box(col="gray80")  
axis(1,cex.axis=0.8,col="gray80",tcl=-0.2)
axis(2,cex.axis=0.8,col="gray80",las=2,tcl=-0.2)

plot(decadedf$year,decadedf$ipi,
     cex=decadedf$propcounts*.3, axes=FALSE,
     xlim=c(2002,2014))
box(col="gray80")
axis(4,cex.axis=0.8,col="gray80",las=2,tcl=-0.2)
abline(mod_f5_ipi,col=alpha("black",0.2),lw=2)

plot(decadedf$year,decadedf$freq,
     cex=decadedf$propcounts*.3, axes=FALSE,
     xlim=c(2002,2014),
     xlab="Year",ylab="Frequency (Hz)")
box(col="gray80")
axis(4,cex.axis=0.8,col="gray80",las=2,tcl=-0.2)
axis(1,cex.axis=0.8,col="gray80",tcl=-0.2)
abline(mod_f5_freq,col=alpha("black",0.2),lw=2)

dev.off()
