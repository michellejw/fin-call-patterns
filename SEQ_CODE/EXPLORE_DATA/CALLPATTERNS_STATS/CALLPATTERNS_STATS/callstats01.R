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






#---------- Figure 5: 2 panels, Decadal F/IPI -----------

# Subset dataframe to include only Decadal analysis stations
decadedf = subset(dat,(station=='AX' | station=='KENE' | station=='KEMF') 
                  & peakcounts > 50)
# Further subset to include only notes in the low IPI/high freq group
decadedf_sub = subset(decadedf,(freq<21)&(ipi>22)) 
# Regression through decadal frequency
mod_f5_freq <- lm(freq~year, data=decadedf_sub)
mod_f5_ipi <- lm(ipi~year, data=decadedf_sub)

flims = c(17.5,24)
ipilims = c(10,36)

# Get legend point sizes
pcount_1 = decadedf$propcounts
pcount_2 = decadedf$peakcounts
count_mid = round(min(pcount_2)+(max(pcount_2)-min(pcount_2))/2)
count_min = min(pcount_2)
count_max = max(pcount_2)
pcount_mid = min(pcount_1)+(max(pcount_1)-min(pcount_1))/2
pcount_min = min(pcount_1)
pcount_max = max(pcount_1)

pdf(file="test.pdf",width=10,height=7)
fig5 = layout(matrix(c(1,2,1,3),2,2,byrow = TRUE))
par(oma=c(5,2,2,2), mgp=c(1.4,.4,0), cex.axis=1.2, cex.lab=1.2)

# Panel 1
par(mar=c(4,4,4,0.2))
plot(decadedf$freq,decadedf$ipi,
     cex=decadedf$propcounts*.6,
     col=cpal_sta[decadedf$netvec],axes=FALSE,
     xlim=flims,ylim=ipilims,
     xlab="",ylab="")
box(col="gray80")  
axis(1,col="gray80",tcl=-0.2)
axis(2,col="gray80",las=2,tcl=-0.2)
mtext("IPI (s)", side=2, line=2)
mtext("Frequency (Hz)", side=1, line=2)
legend(18,14,c(toString(count_min),toString(count_mid),toString(count_max)),
       lwd='NONE',
       bty='n',
       pt.cex=decadedf$propcounts*.6,
       pch=1)


# Panel 2
par(mar=c(0,0.2,4,4))
plot(decadedf$year,decadedf$ipi,
     cex=decadedf$propcounts*.6, axes=FALSE,
     col=cpal_sta[decadedf$netvec],
     xlim=c(2002,2014))
box(col="gray80")
axis(4,col="gray80",las=2,tcl=-0.2)
abline(mod_f5_ipi,col=alpha("black",0.2),lw=2)
mtext("IPI (s)", side=4, line=2)

# Panel 3
par(mar=c(4,0.2,0,4))
plot(decadedf$year,decadedf$freq,
     cex=decadedf$propcounts*.6, axes=FALSE,
     col=cpal_sta[decadedf$netvec],
     xlim=c(2002,2014),
     xlab="",ylab="")
box(col="gray80")
axis(4,col="gray80",las=2,tcl=-0.2)
axis(1,col="gray80",tcl=-0.2)
abline(mod_f5_freq,col=alpha("black",0.2),lw=2)
mtext("Frequency (Hz)", side=4, line=2)
mtext("Year", side=1, line=2)

dev.off()
