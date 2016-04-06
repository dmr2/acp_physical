#!/usr/bin/R

name <- "IPCCAR5climsens_rcp85_DAT_SURFACE_TEMP_BO_15Nov2013_185227.OUT_2080-2099_perc"
diri <- "/Users/dmr/Consulting/RB/magicc"

fil <- paste(diri,name,sep="/")
df <- read.table(fil, header = TRUE, sep = "\t")

par(mar=c(5,6,6,6)+0.1)
plot(df$MAGICC.runs,df$q99,type="n",bty="n",ylab="degrees Celsius above 1981-2010 \n (2080-2099)",xlab="Number of MAGICC runs",ylim=c(0,10), col="blue", cex.lab=1, cex.axis=1, cex.main=1, cex.sub=1,
     main="")
lines(df$MAGICC.runs,df$q99,lwd=3,col=("red"))
lines(df$MAGICC.runs,df$q95,lwd=3,col=("orange"))
lines(df$MAGICC.runs,df$q83,lwd=3,col=("yellow"))
lines(df$MAGICC.runs,df$q50,lwd=3,col=("green"))
lines(df$MAGICC.runs,df$q17,lwd=3,col=("violet"))
lines(df$MAGICC.runs,df$q5,lwd=3,col=("blue"))
