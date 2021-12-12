#!/usr/bin/env Rscript

#
#------------------------------------------------------------------------------
# ~-~ plot_tpwave.r -~-
# Plot the output of tpwave.r
# 
# Author: CL (cristianl@met.no)
#------------------------------------------------------------------------------

library(raster)
library(argparser)

options(warn=2)

#
#------------------------------------------------------------------------------
# Read command line arguments

p <- arg_parser("tpwave")
#..............................................................................
p <- add_argument(p, "--extent_label",
                  help="domain label",
                  type="character",
                  default=NA)
p <- add_argument(p, "--season_label",
                  help="season label (djf,mam,jja,son)",
                  type="character",
                  default=NA)
p <- add_argument(p, "--model_lab",
                  help="model label",
                  type="character",
                  default="era5")
p <- add_argument(p, "--filter_len",
                  help="filter length (days)",
                  type="integer",
                  default=365)
#..............................................................................
argv <- parse_args(p)

#------------------------------------------------------------------------------
# Define the correct file name(s)

if (argv$model_lab == "era5") {
  dir <- "./data/"
  if ( is.na(argv$extent_label)) {
    pattern <- paste0( "ERA5_tp_day_waveEn_ALL*")
  } else {
    pattern <- paste0( "ERA5_tp_day_waveEn_", argv$extent_label)
  }
  files <- list.files( dir, pattern = pattern, full.names=T)
} else if (argv$model_lab == "era20c") {
  files <- "out_era20c.rdata"
} else if (argv$model_lab == "nmc48") {
  files <- "out_nmc48.rdata"
}

#
#------------------------------------------------------------------------------
# Get some dates in R format

Rdate51 <- strptime( "1951-01-01T00", format="%Y-%m-%dT%H", tz="UTC")
Rdate61 <- strptime( "1961-01-01T00", format="%Y-%m-%dT%H", tz="UTC")
Rdate80 <- strptime( "1980-12-31T23", format="%Y-%m-%dT%H", tz="UTC")
Rdate90 <- strptime( "1990-12-31T23", format="%Y-%m-%dT%H", tz="UTC")
Rdate91 <- strptime( "1991-01-01T00", format="%Y-%m-%dT%H", tz="UTC")
Rdate20 <- strptime( "2020-12-31T23", format="%Y-%m-%dT%H", tz="UTC")

#
#------------------------------------------------------------------------------
# Define the right subset

str <- argv$extent_label 
if ( is.na(argv$extent_label)) str <- "ALL"
if ( is.na(argv$season_label)) {
  str <- paste0( str, "_all")
  mm_in <- 1:12
} else {
  str <- paste0( str, "_", argv$season_label) 
  if ( argv$season_label == "djf") mm_in <- c( 12,  1, 2)
  if ( argv$season_label == "mam") mm_in <- c(  3,  4, 5)
  if ( argv$season_label == "jja") mm_in <- c(  6,  7, 8)
  if ( argv$season_label == "son") mm_in <- c(  9, 10,11)
}

#
#------------------------------------------------------------------------------
# Read input data

ts       <- integer(0)
which_mx <- integer(0)
res <- NA
a<-vector()
a_l <- vector()

for (file in files) {
  
  load( file)
  
  # first time in: define data structures
  if (any(is.na(res))) {
    res <- res(s)
    j_fw <- dim(En2o_t)[1]
    j_mw <- dim(En2o_t)[1]-1
    for (l in 1:j_mw) {
      a[l] <- 2**(l-1) * res(s)[1] * 2**(l-1) * res(s)[2]
      a_l[l] <- sqrt( a[l])
      if (a_l[l]<1) { a_l[l] <- round(a_l[l],1) } else { a_l[l] <- round(a_l[l]) }
    }
    a[j_fw] <- 2**(j_fw-1) * res(s)[1] * 2**(j_fw-1) * res(s)[2]
    a_l[j_fw] <- round(sqrt( a[j_fw]))
  }

  # month
  mm <- as.POSIXlt(tseq[ixt],tz="UTC")$mon+1

  # select dates
  ixmm <- which( ((as.POSIXlt(tseq[ixt],tz="UTC")$mon+1) %in% mm_in) & ((as.POSIXlt(tseq[ixt],tz="UTC")$year+1900) >= 1950))
  ix_nona <- integer(0)
  for (i in ixmm) {
    if ( !any( is.na( En2o_t[,i]))) ix_nona <- c( ix_nona, i)
  }
  ixmm <- ix_nona
  
  # compute energies
  En2o_t <- En2o_t[,ixmm]
  En2o_t.tot_all <- colSums(En2o_t)
  En2o_t.tot_mw  <- colSums(En2o_t[1:j_mw,])
 
  En2o_t.perc_all <- En2o_t; En2o_t.perc_all[] <- NA
  En2o_t.perc_mw  <- En2o_t; En2o_t.perc_mw[]  <- NA
  for (i in 1:length(En2o_t.tot_all))
    En2o_t.perc_all[1:j_fw,i] <- 100 * En2o_t[1:j_fw,i] / En2o_t.tot_all[i]
  for (i in 1:length(En2o_t.tot_all)) 
    En2o_t.perc_mw[1:j_mw,i]  <- 100 * En2o_t[1:j_mw,i] / En2o_t.tot_mw[i]
  
  # define data structures used in plotting figures
  if ( !exists( "En2_t")) {
    En2_t          <- En2o_t
    En2_t.perc_all <- En2o_t.perc_all
    En2_t.perc_mw  <- En2o_t.perc_mw
  } else {
    En2_t          <- cbind( En2_t, En2o_t)
    En2_t.perc_all <- cbind( En2_t.perc_all, En2o_t.perc_all)
    En2_t.perc_mw  <- cbind( En2_t.perc_mw,  En2o_t.perc_mw)
  }
  ts    <- c( ts, tseq[ixt[ixmm]])

  # get the spatial level where the energy is max
  aux <- vector( mode="numeric", length=length(ixt[ixmm]))
  for (i in 1:length(ixmm)) {
    aux[i] <- which.max( En2o_t[1:j_mw,i])
  }
  which_mx <- c( which_mx, aux)
} # end loop over files

#
#------------------------------------------------------------------------------
# Fig. Plot all squared energies

nt<- length(ts)
range <- c(0,max(En2_t[1:j_mw,],na.rm=T))

ffout <- paste0( "en2_hist_",str,"_1950-2020.png")

png( file=ffout, width=800, height=800)
plot( 1:j_mw, En2_t[1:j_mw,1], ylim=range, col="white", axes=F, xlab="Spatial scale [deg]", ylab="Squared Energy")
for( i in 1:nt) {
  lines( 1:j_mw, En2_t[1:j_mw,i], col="gray")
}
lines( 1:j_mw, rowMeans(En2_t[1:j_mw,]), lwd=5, col="black")
axis(2)
axis(1,at=1:j_mw,labels=a_l[1:j_mw])
aux <- dev.off()
cat( paste("written file",ffout,"\n"))

#
#------------------------------------------------------------------------------
# Fig. Plot squared energies for 1961-1990

ix1 <- which( as.numeric(ts) >= as.numeric(Rdate61) & 
              as.numeric(ts) <= as.numeric(Rdate90))
ix2 <- which( as.numeric(ts) >= as.numeric(Rdate91) & 
              as.numeric(ts) <= as.numeric(Rdate20))

nt<- length(ix1)
range <- c(0, max(En2_t[1:j_mw,c(ix1,ix2)],na.rm=T))

en1 <- rowMeans(En2_t[1:j_mw,ix1])
en1qs <- apply( En2_t[1:j_mw,ix1], MAR=1, FUN=function(x){as.numeric(quantile(x,probs=0.01))})
en1qb <- apply( En2_t[1:j_mw,ix1], MAR=1, FUN=function(x){as.numeric(quantile(x,probs=0.99))})
en2 <- rowMeans(En2_t[1:j_mw,ix2])
en2qs <- apply( En2_t[1:j_mw,ix2], MAR=1, FUN=function(x){as.numeric(quantile(x,probs=0.01))})
en2qb <- apply( En2_t[1:j_mw,ix2], MAR=1, FUN=function(x){as.numeric(quantile(x,probs=0.99))})

ffout <- paste0( "en2_hist_",str,"_1961-1990.png")
png( file=ffout, width=800, height=800)
plot( 1:j_mw, En2_t[1:j_mw,1], ylim=range, col="white", axes=F,
      xlab="Spatial scale [deg]", ylab="Squared Energy",
      main="1961-1990", cex.main=2)
for( i in 1:nt) {
  lines( 1:j_mw, En2_t[1:j_mw,ix1[i]], col="gray")
}
lines( 1:j_mw, en1, lwd=5, col="black")
lines( 1:j_mw, en1qs, lwd=3, col="black", lty=2)
lines( 1:j_mw, en1qb, lwd=3, col="black", lty=2)
axis(2)
axis(1,at=1:j_mw,labels=a_l[1:j_mw])
aux <- dev.off()
cat( paste("written file",ffout,"\n"))

#
#------------------------------------------------------------------------------
# Fig. Plot squared energies for 1991-2020

nt<- length(ix2)

ffout <- paste0( "en2_hist_",str,"_1991-2020.png")
png( file=ffout, width=800, height=800)
plot( 1:j_mw, En2_t[1:j_mw,1], ylim=range, col="white", axes=F,
      xlab="Spatial scale [deg]", ylab="Squared Energy",
      main="1991-2020", cex.main=2)
for( i in 1:nt) {
  lines( 1:j_mw, En2_t[1:j_mw,ix2[i]], col="gray")
}
lines( 1:j_mw, en2, lwd=5, col="black")
lines( 1:j_mw, en2qs, lwd=3, col="black", lty=2)
lines( 1:j_mw, en2qb, lwd=3, col="black", lty=2)
axis(2)
axis(1,at=1:j_mw,labels=a_l[1:j_mw])
aux <- dev.off()
cat( paste("written file",ffout,"\n"))

#
#------------------------------------------------------------------------------
# Fig. at which decomposition level is the max of the squared energy?

vl <- seq( strptime("1950-01-01",format="%Y-%m-%d"), strptime("2020-01-01",format="%Y-%m-%d"), by="1 year")
at <- seq( strptime("1950-06-01",format="%Y-%m-%d"), strptime("2020-06-01",format="%Y-%m-%d"), by="10 year")

ffout <- paste0( "en2_ts_mode_",str,"_1950-2020.png")
png( file=ffout, width=800, height=800)
par(mar=c(5,5,1,5))
plot( ts, which_mx, col="gray", type="l", axes=F, xlab="", ylab="")
axis(2)
mtext(side=1,line=3,text="Year",cex=2)
mtext(side=2,line=3,text="level, depth of the decomposition (daily, gray)",cex=2)
par(new=T)
par(mar=c(5,5,1,5))
plot( ts, filter( which_mx, rep(1/(argv$filter_len*5),(argv$filter_len*5))), type="l",axes=F,lwd=5, xlab="", ylab="")
abline(v=vl,col="black",lty=2)
axis(1,at=at,labels=format(at,"%Y"))
axis(4)
mtext(side=4,line=3,text="level, depth of the decomposition (5-y mean, black)",cex=2)
aux <- dev.off()
cat( paste("written file",ffout,"\n"))

#
#------------------------------------------------------------------------------
# Fig. squared energy, comparison of the spectrum over different periods

Rdate91 <- strptime( "1991-01-01T00", format="%Y-%m-%dT%H", tz="UTC")
Rdate20 <- strptime( "2020-12-31T23", format="%Y-%m-%dT%H", tz="UTC")
ix1 <- which( as.numeric(ts) >= as.numeric(Rdate51) & 
              as.numeric(ts) <= as.numeric(Rdate80))
ix2 <- which( as.numeric(ts) >= as.numeric(Rdate91) & 
              as.numeric(ts) <= as.numeric(Rdate20))
nt1<- length(ts)
en1 <- rowMeans(En2_t[1:j_mw,ix1])
en1qs <- apply( En2_t[1:j_mw,ix1], MAR=1, FUN=function(x){as.numeric(quantile(x,probs=0.01))})
en1qb <- apply( En2_t[1:j_mw,ix1], MAR=1, FUN=function(x){as.numeric(quantile(x,probs=0.99))})
en2 <- rowMeans(En2_t[1:j_mw,ix2])
en2qs <- apply( En2_t[1:j_mw,ix2], MAR=1, FUN=function(x){as.numeric(quantile(x,probs=0.01))})
en2qb <- apply( En2_t[1:j_mw,ix2], MAR=1, FUN=function(x){as.numeric(quantile(x,probs=0.99))})
range <- c(0,max(c(en1qb,en2qb),na.rm=T))
ffout <- paste0( "en2_histcomp_",str,"_1951-1980vs1991-2020.png")
png( file=ffout, width=800, height=800)
par(mar=c(5,5,1,1))
plot( 1:j_mw, En2_t[1:j_mw,1], ylim=range, col="white", axes=F, xlab="", ylab="")
polygon(c(1:j_mw,j_mw:1), c(en1qs,en1qb[j_mw:1]),density=10,angle=45,col="cornflowerblue")
polygon(c(1:j_mw,j_mw:1), c(en2qs,en2qb[j_mw:1]),density=10,col="red",angle=-45)
lines( 1:j_mw, en1, lwd=12, col="blue")
lines( 1:j_mw, en2, lwd=8, col="red")
axis(2,cex.axis=2)
axis(1,at=1:j_mw,labels=a_l[1:j_mw],cex.axis=2)
abline(v=1:j_mw,lty=2,col="gray")
abline(h=0,lty=1)
mtext(side=1,line=3,text="spatial scale (degrees)",cex=2)
mtext(side=2,line=3,text="squared energy (mm^2)",cex=2)
aux <- dev.off()
cat( paste("written file",ffout,"\n"))

ix1 <- which( as.numeric(ts) >= as.numeric(Rdate61) & 
              as.numeric(ts) <= as.numeric(Rdate90))
ix2 <- which( as.numeric(ts) >= as.numeric(Rdate91) & 
              as.numeric(ts) <= as.numeric(Rdate20))
nt1<- length(ts)
en1 <- rowMeans(En2_t[1:j_mw,ix1])
en1qs <- apply( En2_t[1:j_mw,ix1], MAR=1, FUN=function(x){as.numeric(quantile(x,probs=0.01))})
en1qb <- apply( En2_t[1:j_mw,ix1], MAR=1, FUN=function(x){as.numeric(quantile(x,probs=0.99))})
en2 <- rowMeans(En2_t[1:j_mw,ix2])
en2qs <- apply( En2_t[1:j_mw,ix2], MAR=1, FUN=function(x){as.numeric(quantile(x,probs=0.01))})
en2qb <- apply( En2_t[1:j_mw,ix2], MAR=1, FUN=function(x){as.numeric(quantile(x,probs=0.99))})
range <- c(0,max(c(en1qb,en2qb),na.rm=T))
range <- c(0,max(En2_t,na.rm=T))
ffout <- paste0( "en2_histcomp_",str,"_1961-1990vs1991-2020.png")
png( file=ffout, width=800, height=800)
par(mar=c(5,5,1,1))
plot( 1:j_mw, En2_t[1:j_mw,1], ylim=range, col="white", axes=F, xlab="", ylab="")
polygon(c(1:j_mw,j_mw:1), c(en1qs,en1qb[j_mw:1]),density=10,angle=45,col="cornflowerblue")
polygon(c(1:j_mw,j_mw:1), c(en2qs,en2qb[j_mw:1]),density=10,col="red",angle=-45)
#boxplot( 1:j_mw, En2_t[1:j_mw,1])
lines( 1:j_mw, en1, lwd=12, col="blue")
lines( 1:j_mw, en2, lwd=8, col="red")
axis(2,cex.axis=2)
axis(1,at=1:j_mw,labels=a_l[1:j_mw],cex.axis=2)
abline(v=1:j_mw,lty=2,col="gray")
abline(h=0,lty=1)
mtext(side=1,line=3,text="spatial scale (degrees)",cex=2)
mtext(side=2,line=3,text="squared energy (mm^2)",cex=2)
aux <- dev.off()
cat( paste("written file",ffout,"\n"))

#
#------------------------------------------------------------------------------
# Fig. squared energy, comparison of the spectrum over different periods (boxplots)
# (figure in the paper)

range <- c(0,max(En2_t[,c(ix1,ix2)],na.rm=T))
ffout <- paste0( "en2_bxp_",str,"_1961-1990vs1991-2020.png")
png( file=ffout, width=800, height=800)
par(mar=c(5,5,1,1))
plot( 1:j_mw, En2_t[1:j_mw,1], xlim=c(0,j_mw), ylim=range, col="white", axes=F, xlab="", ylab="")
for (i in 1:j_mw) lines(c(i-0.75,i-0.75),c(en1qs[i],en1qb[i]),lty=2,col="cornflowerblue",lwd=3)
for (i in 1:j_mw) points(c(i-0.75,i-0.75),c(en1qs[i],en1qb[i]),pch=3,col="cornflowerblue",cex=2,lwd=3)
for (i in 1:j_mw) lines(c(i-0.25,i-0.25),c(en2qs[i],en2qb[i]),lty=2,col="pink",lwd=3)
for (i in 1:j_mw) points(c(i-0.25,i-0.25),c(en2qs[i],en2qb[i]),pch=3,col="pink",cex=2,lwd=3)
for (i in 1:j_mw) boxplot( En2_t[i,ix2], add=T, at=i-0.25, col="pink", range=0.000001, outline=F, axes=F, boxwex=0.6)
for (i in 1:j_mw) boxplot( En2_t[i,ix1], add=T, at=i-0.75, col="cornflowerblue", range=0.000001, outline=F, axes=F, boxwex=0.6)
for (i in 1:j_mw) {
  ix_out <- which( En2_t[i,ix2] > en2qb[i] | En2_t[i,ix2] < en2qs[i])
  if ( length( ix_out)>0) points( rep( i-0.25, length( ix_out)), En2_t[i,ix2[ix_out]],cex=1,col="pink") 
  ix_out <- which( En2_t[i,ix1] > en1qb[i] | En2_t[i,ix1] < en1qs[i])
  if ( length( ix_out)>0) points( rep( i-0.75, length( ix_out)), En2_t[i,ix1[ix_out]],cex=1,col="cornflowerblue") 
}
abline(v=1:j_mw,lty=2,col="gray")
axis(2,cex.axis=2)
axis(1,at=((1:j_mw)-0.5),labels=a_l[1:j_mw],cex.axis=2)
abline(h=0,lty=1)
mtext(side=1,line=3,text="spatial scale (degrees)",cex=2)
mtext(side=2,line=3,text="squared energy (mm^2)",cex=2)
range <- c(0,max(c(en1qb,en2qb),na.rm=T))
par(new=T)
par(mar=c(38,38,1,1))
plot( 1:j_mw, En2_t[1:j_mw,1], ylim=range, col="white", axes=F, xlab="", ylab="")
rect(-j_mw,-100,2*j_mw,1000,col="white")
polygon(c(1:j_mw,j_mw:1), c(en1qs,en1qb[j_mw:1]),density=10,angle=45,col="cornflowerblue")
polygon(c(1:j_mw,j_mw:1), c(en2qs,en2qb[j_mw:1]),density=10,col="red",angle=-45)
#boxplot( 1:j_mw, En2_t[1:j_mw,1])
lines( 1:j_mw, en1, lwd=12, col="blue")
lines( 1:j_mw, en2, lwd=8, col="red")
abline(h=0,v=1)
box()
aux <- dev.off()

#
#------------------------------------------------------------------------------
# Fig. squared energy as percentage of the total energy, comparison of the 
#       spectrum over different periods

en1 <- rowMeans(En2_t.perc_all[1:j_mw,ix1])
en1qs <- apply( En2_t.perc_all[1:j_mw,ix1], MAR=1, FUN=function(x){as.numeric(quantile(x,probs=0.01))})
en1qb <- apply( En2_t.perc_all[1:j_mw,ix1], MAR=1, FUN=function(x){as.numeric(quantile(x,probs=0.99))})
en2 <- rowMeans(En2_t.perc_all[1:j_mw,ix2])
en2qs <- apply( En2_t.perc_all[1:j_mw,ix2], MAR=1, FUN=function(x){as.numeric(quantile(x,probs=0.01))})
en2qb <- apply( En2_t.perc_all[1:j_mw,ix2], MAR=1, FUN=function(x){as.numeric(quantile(x,probs=0.99))})
range <- c( 0, max(en1qb,en2qb,na.rm=T))
ffout <- paste0( "en2pall_histcomp_",str,"_1961-1990vs1991-2020.png")
png( file=ffout, width=800, height=800)
par(mar=c(5,5,1,1))
plot( 1:j_mw, En2_t.perc_all[1:j_mw,1], ylim=range, col="white", axes=F, xlab="", ylab="")
polygon(c(1:j_mw,j_mw:1), c(en1qs,en1qb[j_mw:1]),density=10,angle=45,col="cornflowerblue")
polygon(c(1:j_mw,j_mw:1), c(en2qs,en2qb[j_mw:1]),density=10,col="red",angle=-45)
lines( 1:j_mw, en1, lwd=12, col="blue")
lines( 1:j_mw, en2, lwd=8, col="red")
axis(2,cex.axis=2)
axis(1,at=1:j_mw,labels=a_l[1:j_mw],cex.axis=2)
abline(v=1:j_mw,lty=2,col="gray")
abline(h=0,lty=1)
mtext(side=1,line=3,text="spatial scale (degrees)",cex=2)
mtext(side=2,line=3,text="percentage of the total(+FW) squared energy (%)",cex=2)
aux <- dev.off()
cat( paste("written file",ffout,"\n"))

en1 <- rowMeans(En2_t.perc_mw[1:j_mw,ix1])
en1qs <- apply( En2_t.perc_mw[1:j_mw,ix1], MAR=1, FUN=function(x){as.numeric(quantile(x,probs=0.01))})
en1qb <- apply( En2_t.perc_mw[1:j_mw,ix1], MAR=1, FUN=function(x){as.numeric(quantile(x,probs=0.99))})
en2 <- rowMeans(En2_t.perc_mw[1:j_mw,ix2])
en2qs <- apply( En2_t.perc_mw[1:j_mw,ix2], MAR=1, FUN=function(x){as.numeric(quantile(x,probs=0.01))})
en2qb <- apply( En2_t.perc_mw[1:j_mw,ix2], MAR=1, FUN=function(x){as.numeric(quantile(x,probs=0.99))})
range <- c( 0, max(en1qb,en2qb,na.rm=T))
#range <- range(c(en1qb,en2qb),na.rm=T)
ffout <- paste0( "en2pmw_histcomp_",str,"_1961-1990vs1991-2020.png")
png( file=ffout, width=800, height=800)
par(mar=c(5,5,1,1))
plot( 1:j_mw, En2_t.perc_mw[1:j_mw,1], ylim=range, col="white", axes=F, xlab="", ylab="")
polygon(c(1:j_mw,j_mw:1), c(en1qs,en1qb[j_mw:1]),density=10,angle=45,col="cornflowerblue")
polygon(c(1:j_mw,j_mw:1), c(en2qs,en2qb[j_mw:1]),density=10,col="red",angle=-45)
lines( 1:j_mw, en1, lwd=12, col="blue")
lines( 1:j_mw, en2, lwd=8, col="red")
axis(2,cex.axis=2)
axis(1,at=1:j_mw,labels=a_l[1:j_mw],cex.axis=2)
abline(v=1:j_mw,lty=2,col="gray")
abline(h=0,lty=1)
mtext(side=1,line=3,text="spatial scale (degrees)",cex=2)
mtext(side=2,line=3,text="percentage of the total squared energy (%)",cex=2)
aux <- dev.off()
cat( paste("written file",ffout,"\n"))

#
#------------------------------------------------------------------------------
# Fig. time series of the squared energies (level by level)

at <- seq( strptime("1950-06-01",format="%Y-%m-%d"), strptime("2020-06-01",format="%Y-%m-%d"), by="10 year")
for (j in 1:j_fw) {
  jstr <- paste0("MW",formatC(j,flag="0",width=2)); if (j==j_fw) jstr <- "FW" 
  ffout <- paste0( "tpmean2_ts_",str,"_",jstr,"_1950-2020.png")
  png( file=ffout, width=800, height=800)
  par( mar=c(5,5,1,1))
  plot( ts, filter( En2_t[j,], rep(1/(5*argv$filter_len),(5*argv$filter_len))), type="l",lwd=5, axes=F, xlab="", ylab="")
  mtext(side=1,line=3,text="Year",cex=2)
  mtext(side=2,line=3,"precipitation energy (mm^2) (5-y mean)", cex=2)
  mtext(side=3,line=-1.3,text=paste("j=",jstr," L=",a_l[j]," deg"),cex=3)
  axis(2,cex.axis=2)
  axis(1,at=at,labels=format(at,"%Y"),cex.axis=2)
  aux <- dev.off()
  cat( paste("written file",ffout,"\n"))
}

#
#------------------------------------------------------------------------------
# Fig. time series of the squared energies % of all energy (level by level)

at <- seq( strptime("1950-06-01",format="%Y-%m-%d"), strptime("2020-06-01",format="%Y-%m-%d"), by="10 year")
for (j in 1:j_fw) {
  jstr <- paste0("MW",formatC(j,flag="0",width=2)); if (j==j_fw) jstr <- "FW" 
  ffout <- paste0( "tpmean2pall_ts_",str,"_",jstr,"_1950-2020.png")
  png( file=ffout, width=800, height=800)
  par( mar=c(5,5,1,1))
  plot( ts, filter( En2_t.perc_all[j,], rep(1/(5*argv$filter_len),(5*argv$filter_len))), type="l",lwd=5, axes=F, xlab="", ylab="")
  mtext(side=1,line=3,text="Year",cex=2)
  mtext(side=2,line=3,"percentage of the total(+FW) squared energy (%) (5-y mean)", cex=2)
  mtext(side=3,line=-1.3,text=paste("j=",jstr," L=",a_l[j]," deg"),cex=3)
  axis(2,cex.axis=2)
  axis(1,at=at,labels=format(at,"%Y"),cex.axis=2)
  aux <- dev.off()
  cat( paste("written file",ffout,"\n"))
}

#
#------------------------------------------------------------------------------
# Fig. time series of the squared energies % of total energy in wavelets (level by level)

at <- seq( strptime("1950-06-01",format="%Y-%m-%d"), strptime("2020-06-01",format="%Y-%m-%d"), by="10 year")
for (j in 1:j_mw) {
  jstr <- paste0("MW",formatC(j,flag="0",width=2)); if (j==j_fw) jstr <- "FW" 
  ffout <- paste0( "tpmean2pmw_ts_",str,"_",jstr,"_1950-2020.png")
  png( file=ffout, width=800, height=800)
  par( mar=c(5,5,1,1))
  plot( ts, filter( En2_t.perc_mw[j,], rep(1/(5*argv$filter_len),(5*argv$filter_len))), type="l",lwd=5, axes=F, xlab="", ylab="")
  mtext(side=1,line=3,text="Year",cex=2)
  mtext(side=2,line=3,"percentage of the total squared energy (%) (5-y mean)", cex=2)
  mtext(side=3,line=-1.3,text=paste("j=",jstr," L=",a_l[j]," deg"),cex=3)
  axis(2,cex.axis=2)
  axis(1,at=at,labels=format(at,"%Y"),cex.axis=2)
  aux <- dev.off()
  cat( paste("written file",ffout,"\n"))
}
