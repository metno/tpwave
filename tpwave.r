#!/usr/bin/env Rscript

library(ncdf4)
library(waveslim)
library(raster)
library(argparser)

#
#------------------------------------------------------------------------------
# Read command line arguments

p <- arg_parser("wave_it")

p <- add_argument(p, "--date1", help="begin (%Y-%m-%d)", type="character", default="1950-01-01")

p <- add_argument(p, "--date2", help="end (%Y-%m-%d)", type="character", default="1950-01-02")

p <- add_argument(p, "--ffin", help="input file", type="character", default="./data/ERA5_tp_day_1950-01-01_1950-01-02.nc")

p <- add_argument(p, "--ffout", help="output file", type="character", default="ERA5_tp_day_example_waveEn.RData")

p <- add_argument(p, "--extent", help="extent", type="numeric", nargs=Inf, default=NA)

p <- add_argument(p, "--extent_lab", help="label", type="character", default=NA)

argv <- parse_args(p)

#
#------------------------------------------------------------------------------
# Define extent based on the command line argurment
if ( !is.na( argv$extent_lab)) {
  if ( argv$extent_lab == "EU") argv$extent <- c( -40,  +75,  10,  85)
  if ( argv$extent_lab == "AF") argv$extent <- c( -40,  +75, -55,  45)
  if ( argv$extent_lab == "AS") argv$extent <- c( +25, +179,  -5,  85)
  if ( argv$extent_lab == "OC") argv$extent <- c( +40, +179, -70,  10)
  if ( argv$extent_lab == "NA") argv$extent <- c(-179,  -10,  10,  85)
  if ( argv$extent_lab == "SA") argv$extent <- c(-179,  -10, -70,  20)
  # nothern hemisphere
  if ( argv$extent_lab == "NH") argv$extent <- c(-179, +179,   0,  85)
  # southern hemisphere
  if ( argv$extent_lab == "SH") argv$extent <- c(-179, +179, -85,   0)
  # tropics
  if ( argv$extent_lab == "TR") argv$extent <- c(-179, +179, -25,  25)
  # temperate north
  if ( argv$extent_lab == "TN") argv$extent <- c(-179, +179,  23,  67)
  # temperate south
  if ( argv$extent_lab == "TS") argv$extent <- c(-179, +179, -67, -23)
}

#
#------------------------------------------------------------------------------
# Constants

date_format <- "%Y-%m-%d"
proj4.wgs84     <- "+proj=longlat +datum=WGS84"

#
#------------------------------------------------------------------------------
# time 

if ( !is.na(argv$date1) & !is.na(argv$date2)) {
  Rdate1 <- strptime( paste0(argv$date1, "T00"), format=paste0( date_format, "T%H"), tz="UTC")
  Rdate2 <- strptime( paste0(argv$date2, "T23"), format=paste0( date_format, "T%H"), tz="UTC")
}

#
#------------------------------------------------------------------------------
# open input file

if ( !file.exists( argv$ffin)) q()
nc <-nc_open( argv$ffin)
tp <-nc$var[[2]]
ndims <- tp$ndims
start <- rep(1,ndims)
varsize <- tp$varsize

#
#------------------------------------------------------------------------------
# define timesteps to read from file

tseq <- as.POSIXct( (tp$dim[[3]]$vals * 60 * 60), origin="1900-01-01", tz="UTC")

if ( !is.na(argv$date1) & !is.na(argv$date2)) {
  ixt <- which( as.numeric(tseq) >= as.numeric(Rdate1) & 
                as.numeric(tseq) <= as.numeric(Rdate2))
} else {
  ixt <- 1:length(tseq)
}

#
#------------------------------------------------------------------------------
# define spatial grids

lons_tot <- tp$dim[[1]]$vals
lats_tot <- tp$dim[[2]]$vals
w <- raster( nrows=length(lats_tot), ncols=2*length(lons_tot),  
             xmn=min(lons_tot-360)-0.125, xmx=max(lons_tot)+0.125, 
             ymn=min(lats_tot)-0.125, ymx=max(lats_tot)+0.125, 
             crs=proj4.wgs84)
r <- crop(w, extent(c(-180.1, 179.8, -90.125, 90.125)))

if ( !any( is.na( argv$extent))) {
  q <- crop( r, argv$extent)
} else {
  q <- r
}

xy<-xyFromCell(q,1:ncell(q))
lons<-sort(unique(xy[,1]))
lats<-sort(unique(xy[,2]),decreasing=T)
rm(xy)

# number of rows (cols) of the dyadic grid 
dimdy  <- 2**ceiling( log2( max( dim(q)[1:2])))

# dyadic grid
s <- raster( nrows=dimdy, ncols=dimdy,  
             xmn=min(lons)-0.125, xmx=max(lons)+0.125,
             ymn=min(lats)-0.125, ymx=max(lats)+0.125, 
             crs=proj4.wgs84)

nnscales <- log2(dimdy)+1
N <- nnscales-1
listscales <- 2**(seq(1,nnscales)-1)*1

if ( length(ixt) == 0) q()

#
#------------------------------------------------------------------------------
# elaboration

nt      <- length(ixt)
En2o_t<-array(data=NA,dim=c(nnscales,nt))

for( i in 1:nt ) {
  if ( ( i %% 10) == 0) cat( paste(i,"/",nt,"\n"))
  start <- rep(1,ndims)
  start[ndims] <- ixt[i]
  count <- varsize
  count[ndims] <- 1
  data <- ncvar_get( nc, tp, start=start, count=count )
  w[] <- cbind( t(data), t(data)) * 1000
  r <- crop(w, extent(c(-180.1, 179.8, -90.125, 90.125)))
  if ( !any( is.na( argv$extent))) {
    q <- crop( r, argv$extent)
  } else {
    q <- r
  }
  u    <- resample(q,s,method="bilinear")
  obs  <- as.matrix(u)
  Eo.dwt<-dwt.2d(obs, wf = "haar", J = N)
  En2o<-vector(mode="numeric",length=N)
  for (j in 1:N) {
    En2o[j] <- mean((Eo.dwt[[1 + 3 * (j - 1)]]/2^j)^2) + 
               mean((Eo.dwt[[2 + 3 * (j - 1)]]/2^j)^2) +
               mean((Eo.dwt[[3 + 3 * (j - 1)]]/2^j)^2)
  }
  En2o_t[1:N,i] <- En2o[1:N]
  En2o_t[(N+1),i]<-cellStats(u,stat="mean",na.rm=T)**2
}
nc_close(nc)

#
#------------------------------------------------------------------------------
# write output file (Rdata-format)
# r= raster with the specification of the original ERA5 grid
# q= raster with the specification of the original ERA5 grid, but cropped
# s= raster with the specification of the high-res dyadic domain
# lons/lats = longitudes and latitudes of the "q" grid
# dimdy = number of rows/columns of the dyadic domain "s"
# N = number of spatial scales used in the decompositions (mother wavelets)
# listscales = spatial scale size (in number of gridpoints)
# nnscales = N+1, total number of spatial scales (mother + father)
# En2o_t = matrix with the results of the decomposition. The dimensions are: nnsclaes x number_of_time_steps. For each time step, the first element of the vector is the swavelet coefficient for the smallest scale, the last is the coefficient of the scaling function
# tseq = sequence of date/time steps
# ixt = indexes of the time steps actually used
#
save( file=argv$ffout, r, q, s, lons, lats, 
                       dimdy, N, listscales, nnscales,
                       En2o_t, tseq, ixt)

#
#------------------------------------------------------------------------------

q()
