#!/usr/bin/env Rscript

library(raster)
library(argparser)

options(warn=2)

#
#------------------------------------------------------------------------------
# Read command line arguments

p <- arg_parser("tpwave")

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

argv <- parse_args(p)

#
#------------------------------------------------------------------------------
# compose the filenames of the files to read

if (argv$model_lab == "era5") {
  dir <- "/home/cristianl/data/ERA5"
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
# Define dates of the two normal periods considered

Rdate61 <- strptime( "1961-01-01T00", format="%Y-%m-%dT%H", tz="UTC")
Rdate90 <- strptime( "1990-12-31T23", format="%Y-%m-%dT%H", tz="UTC")
Rdate91 <- strptime( "1991-01-01T00", format="%Y-%m-%dT%H", tz="UTC")
Rdate20 <- strptime( "2020-12-31T23", format="%Y-%m-%dT%H", tz="UTC")

#
#------------------------------------------------------------------------------
# define variables related to command line arguments

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
# Read data from input files (output of the wavelet decomposition)

ts       <- integer(0)
which_mx <- integer(0)
res <- NA
a   <- vector()
a_l <- vector()

for (file in files) {
  # load input data
  load( file)
  # define spatial scales labels (a_l)
  if (any(is.na(res))) {
    res <- res(s)
    j_fw <- dim(En2o_t)[1]
    j_mw <- dim(En2o_t)[1]-1
    for (l in 1:j_mw) {
      a[l] <- 2**(l-1) * res(s)[1] * 2**(l-1) * res(s)[2]
      a_l[l] <- sqrt( a[l])
      if (a_l[l]<1) { a_l[l] <- round(a_l[l],2) } else { a_l[l] <- round(a_l[l]) }
    }
    a[j_fw] <- 2**(j_fw-1) * res(s)[1] * 2**(j_fw-1) * res(s)[2]
    a_l[j_fw] <- round( sqrt( a[j_fw]))
  }
  # read energies
  mm <- as.POSIXlt(tseq[ixt],tz="UTC")$mon+1
  ixmm <- which( ((as.POSIXlt(tseq[ixt],tz="UTC")$mon+1) %in% mm_in) & ((as.POSIXlt(tseq[ixt],tz="UTC")$year+1900) >= 1950))
  ix_nona <- integer(0)
  for (i in ixmm) {
    if ( !any( is.na( En2o_t[,i]))) ix_nona <- c( ix_nona, i)
  }
  ixmm <- ix_nona
  En2o_t <- En2o_t[,ixmm]
  En2o_t.tot_all <- colSums(En2o_t)
  En2o_t.tot_mw  <- colSums(En2o_t[1:j_mw,])
  # define and fill the definitive data structures with the energies
  En2o_t.perc_all <- En2o_t; En2o_t.perc_all[] <- NA
  En2o_t.perc_mw  <- En2o_t; En2o_t.perc_mw[]  <- NA
  for (i in 1:length(En2o_t.tot_all)) En2o_t.perc_all[1:j_fw,i] <- 100 * En2o_t[1:j_fw,i] / En2o_t.tot_all[i]
  for (i in 1:length(En2o_t.tot_all)) En2o_t.perc_mw[1:j_mw,i]  <- 100 * En2o_t[1:j_mw,i] / En2o_t.tot_mw[i]
  if ( !exists( "En2_t")) {
    En2_t          <- En2o_t
    En2_t.perc_all <- En2o_t.perc_all
    En2_t.perc_mw  <- En2o_t.perc_mw
  } else {
    En2_t          <- cbind( En2_t, En2o_t)
    En2_t.perc_all <- cbind( En2_t.perc_all, En2o_t.perc_all)
    En2_t.perc_mw  <- cbind( En2_t.perc_mw,  En2o_t.perc_mw)
  }
  # time sequence
  ts    <- c( ts, tseq[ixt[ixmm]])
  aux <- vector( mode="numeric", length=length(ixt[ixmm]))
  for (i in 1:length(ixmm)) {
    aux[i] <- which.max( En2o_t[1:j_mw,i])
  }
  which_mx <- c( which_mx, aux)
} # end loop over files

#
#------------------------------------------------------------------------------
# Select the timesteps belonging to the two normal time periods

ix1 <- which( as.numeric(ts) >= as.numeric(Rdate61) & 
              as.numeric(ts) <= as.numeric(Rdate90))
ix2 <- which( as.numeric(ts) >= as.numeric(Rdate91) & 
              as.numeric(ts) <= as.numeric(Rdate20))

nt<- length(ix1)
range <- c(0, max(En2_t[1:j_mw,c(ix1,ix2)],na.rm=T))

#
#------------------------------------------------------------------------------
# We define thre time periods
#  period 1 = 1961--1990 normal time period 
#  period 2 = 1991--2020 normal time period 
#  period 3 = 2021--2050 normal time period 

bxp1 <- list()
bxp2 <- list()
bxp3 <- list()
bxp1_p <- list()
bxp2_p <- list()

mres_reldev_q50 <- vector()
mres_reldev_q50_a <- vector()
mres_reldev_p_q50 <- vector()
mres_reldev_q99 <- vector()
mres_reldev_future_q50 <- vector()
mres_reldev_future_q50_a <- vector()
mres_reldev_future_q99 <- vector()
# loop over the mother wavelet components
for (i in 1:j_mw) {
  # compute percentiles of energy over period 1 considering the 30-years of daily energies  
  bxp1[[i]] <- boxplot( En2_t[i,ix1], plot=F)
  # compute percentiles of energy over period 2 considering the 30-years of daily energies  
  bxp2[[i]] <- boxplot( En2_t[i,ix2], plot=F)
  # add 99-th percentile
  bxp1[[i]]$q99 <- as.numeric( quantile(En2_t[i,ix1], probs=0.99))
  bxp2[[i]]$q99 <- as.numeric( quantile(En2_t[i,ix2], probs=0.99))
  # relative deviation between the medians period2 wrt period1
  mres_reldev_q50[i]   <- bxp2[[i]]$stats[3,1] / bxp1[[i]]$stats[3,1]
  # relative deviation between the medians period1 wrt period2
  mres_reldev_q50_a[i] <- bxp1[[i]]$stats[3,1] / bxp2[[i]]$stats[3,1]
  # relative deviation between the 99-th percentiles period2 wrt period1
  mres_reldev_q99[i]   <- bxp2[[i]]$q99 / bxp1[[i]]$q99
  # compute percentiles of energy (as percentage over the total energy) over period 1
  bxp1_p[[i]] <- boxplot( En2_t.perc_mw[i,ix1], plot=F)
  # compute percentiles of energy (as percentage over the total energy) over period 2
  bxp2_p[[i]] <- boxplot( En2_t.perc_mw[i,ix2], plot=F)
  # relative deviation between the medians period2 wrt period1
  mres_reldev_p_q50[i] <- bxp2_p[[i]]$stats[3,1] / bxp1_p[[i]]$stats[3,1]

  # linear regression coefficients over time period 2
  y <- En2_t[i,ix2]
  lm <- lm(y~ix2)$coefficients

  # simulate energies over time period 3, given the "lm" regression coefficients 
  ix3 <- max(ix2) + 1:(365*30)
  En2_future <- ix3*lm[2]+lm[1] 

  # compute percentiles of energy over period 2 considering the 30-years of daily energies  
  bxp3[[i]] <- boxplot( En2_future, plot=F)
  # relative deviation between the medians period3 wrt period1
  mres_reldev_future_q50[i]   <- bxp3[[i]]$stats[3,1] / bxp1[[i]]$stats[3,1]
  # relative deviation between the medians period3 wrt period2
  mres_reldev_future_q50_a[i] <- bxp3[[i]]$stats[3,1] / bxp2[[i]]$stats[3,1]
  
  # print output (just a selection, it can be easily adapted)
  print( paste( formatC( round( a_l[i],2), width=5,flag="0"), # spatial scale 
                formatC( round( bxp1[[i]]$stats[3,1],3), width=6,flag="0"), # relative deviation between the medians period2 wrt period1
                formatC( round( bxp2[[i]]$stats[3,1],3), width=6,flag="0"), # relative deviation between the medians period2 wrt period1
                formatC( round( mres_reldev_q50_a[i],3), width=6,flag="0"), # relative deviation between the medians period2 wrt period1
                formatC( round( mres_reldev_future_q50[i],3), width=6,flag="0"), # relative deviation between the medians period3 wrt period1
                formatC( round( mres_reldev_future_q50_a[i],3), width=6,flag="0"))) # relative deviation between the medians period3 wrt period2
} # end loop over the mother wavelet components

# Normal exit
q( status=0)
