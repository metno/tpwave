#!/usr/bin/env Rscript

library(ncdf4)
library(waveslim)
library(raster)
library(argparser)
library(rgdal)

#==============================================================================
# FUNCTIONS - FUNCTIONS - FUNCTIONS - FUNCTIONS - FUNCTIONS - FUNCTIONS -
#==============================================================================

#+
plot_multires <- function( env, prefix="fig", shp=NA, par=NA) {
#------------------------------------------------------------------------------
  options( warn = 2)

  if ( !any(is.na(shp))) {
    b_or <- readOGR( shp[1], shp[2], verbose=F)
    b    <- spTransform( b_or, crs(env$r))
  }

  # set dyadic domain
  xmn <- extent( env$r)[1]
  xmx <- extent( env$r)[2]
  ymn <- extent( env$r)[3]
  ymx <- extent( env$r)[4]

  nx <- ncol( env$r)
  ny <- nrow( env$r)

  n <- ceiling( log( max(nx,ny), 2))

  rdyad <- raster( extent( xmn, xmx, ymn, ymx), 
                   ncol=2**n, nrow=2**n, crs=crs( env))
  rdyad[] <- 0

  cat( paste( "dyadic domain, nx ny dx dy >", ncol(rdyad), nrow(rdyad), round( res(rdyad)[1]), round( res(rdyad)[2]),"\n"))

  dw<-list()
  dw$wf <- "haar" 
  dw$n_levs_mx <- 11
  dw$boundary <- "periodic"

  # Initialization of the wavelet structures
  #  dwt_out. dwt.2d class.
  dwt_out <- dwt.2d( as.matrix(rdyad), wf=dw$wf, J=dw$n_levs_mx, boundary=dw$boundary)
  for (l in 1:dw$n_levs_mx) {
    dwt_out[[3*(l-1)+1]][] <- l*10 + 1 # LHi - wavelet coefficients
    dwt_out[[3*(l-1)+2]][] <- l*10 + 2 # HLi - wavelet coefficients
    dwt_out[[3*(l-1)+3]][] <- l*10 + 3 # HHi - wavelet coefficients
  }

  # scaling function coefficients (father wavelet)
  dwt_out[[3*(dw$n_levs_mx-1)+4]][] <- dw$n_levs_mx*10 + 4 # LL dw$n_levs_mx
  for (i in 1:length(dwt_out)) dwt_out[[i]][] <- 0

  # -- wave -- 
  cat( paste( "wavelet decomposition","\n"))

  # interpolate onto the dyadic grid
  cat( paste( "resample...","\r"))
  s <- env$r
  cat( paste( "resample... ok","\n"))
  s[s<0] <- 0

  # define parameters used for plotting the figure
  if ( any( is.na(par))) {
    rng_or <- range( c(getValues(env$r), getValues(s)))
    br_or  <- c(0,0.1,1,2,5,10,100,280)
    col_or <- c("beige", colorRampPalette(c("azure", "cyan", "cornflowerblue", "blue", "purple"))(length(br_or)-2))
  } else {
    rng_or <- par$rng_or
    br_or  <- par$br_or
    col_or <- par$col_or
  }
  dx <- res(rdyad)[1]*100
  dy <- res(rdyad)[2]*100

  a   <- vector()
  a_l <- vector()
  # transformation operator
  cat( paste( "dwt...","\r"))
  dwt  <- dwt.2d( as.matrix(s), wf=dw$wf, J=dw$n_levs_mx, boundary=dw$boundary)
  cat( paste( "dwt... ok","\n"))

  # what happens when we enflate the energy?
  # coefficients are those used for plotting Fig.8 of Benestad et al. (2022) referring to Hurricane Katrina
  #  they have been obtained using the code "tpwave_multires_ampfact.r"
  #  coeff1 is the ratio between column3 and column2 of Tab.1 
  #  coeff2 is the square root of "mres_reldev_future_q50_a" in "tpwave_multires_ampfact.r"
  # coeff1, used to rescale Katrina's wavelet coefficients over the climate 1961-1991
  coeff1 <- sqrt( c( 0.7700456, 0.7512890, 0.7406982, 0.7515059, 0.7958232, 0.8367114, 0.8771996, 0.9101677, 0.9300676, 0.9820044, 0.9169702))
  # coeff2, used to rescale Katrina's wavelet coefficients over the climate 2021-2050
  coeff2 <- sqrt( c( 1.290512, 1.329671, 1.343901, 1.346885, 1.306342, 1.253900, 1.209975, 1.165373, 1.193174, 1.136975, 1.214505))

  dwt_out1 <- dwt_out
  dwt_out2 <- dwt_out

  for (i in 1:length(dwt_out1)) dwt_out1[[i]][] <- 0
  for (i in 1:length(dwt_out2)) dwt_out2[[i]][] <- 0

  # rescale wavelet coefficients
  for (l in 1:dw$n_levs_mx) {
    dwt_out1[[1+3*(l-1)]] <- coeff1[l] * dwt[[1+3*(l-1)]]
    dwt_out1[[2+3*(l-1)]] <- coeff1[l] * dwt[[2+3*(l-1)]]
    dwt_out1[[3+3*(l-1)]] <- coeff1[l] * dwt[[3+3*(l-1)]]
    dwt_out2[[1+3*(l-1)]] <- coeff2[l] * dwt[[1+3*(l-1)]]
    dwt_out2[[2+3*(l-1)]] <- coeff2[l] * dwt[[2+3*(l-1)]]
    dwt_out2[[3+3*(l-1)]] <- coeff2[l] * dwt[[3+3*(l-1)]]
  }
  
  dwt_out1[[4+3*(dw$n_levs_mx-1)]] <- dwt[[4+3*(dw$n_levs_mx-1)]]
  dwt_out2[[4+3*(dw$n_levs_mx-1)]] <- dwt[[4+3*(dw$n_levs_mx-1)]]

  qqor   <- rdyad
  qqin1   <- rdyad
  qqin2   <- rdyad

#  dwt[[4+3*(dw$n_levs_mx-1)]][] <- 0

  qqor[] <- idwt.2d( dwt)
  qqin1[] <- idwt.2d( dwt_out1)
  qqin2[] <- idwt.2d( dwt_out2)

#  x1 <- -107
#  x2 <- -70
#  y1 <- 10
#  y2 <- 34
  x1 <- argv$extent[1]
  x2 <- argv$extent[2]
  y1 <- argv$extent[3]
  y2 <- argv$extent[4]

  qqor<-crop(qqor,c(x1,x2,y1,y2))
  qqin1<-crop(qqin1,c(x1,x2,y1,y2))
  qqin2<-crop(qqin2,c(x1,x2,y1,y2))
  qqor_a <- qqor
  qqin1_a <- qqin1
  qqin2_a <- qqin2
  minn <- min( getValues(qqor))
  maxx <- max( getValues(qqor))
  maxx <- as.numeric( quantile( getValues(qqor), probs=0.999))
  maxx <- 350
  qqor[qqor>maxx] <- maxx
  qqor[qqor<1] <- 0  
  qqin1[qqin1<1] <- 0  
  qqin2[qqin2<1] <- 0 
  qqin1[qqin1>maxx] <- maxx
  qqin2[qqin2>maxx] <- maxx

  # call the library used to load nice color tables
#  library(fool)
#  col_or<-load_color_table(path="/home/cristianl/projects/rpackages/fool/color_tables",abbrv="precip_11lev")
  # can be replaced with the following
  col_or <- c( "beige",rev( rainbow(11)))
  br_or  <- seq( minn, maxx, length=(length(col_or)+1))
  br_or <- c(0, 1, 5, 10, 20, 30, 40, 50, 70, 120, 180, 240, 300)

  # plot the panel (b) of Fig. 8 in Benestad et al. (2022)
  # original precipitation field
  png( file=paste0( "katrina_orig_2005.png"), height=800, width=1200)
  par(mar=c(3,3,1,1))
  image( qqor, main="", xlab="", ylab="", col=col_or, breaks=br_or,xlim=c(x1,x2),ylim=c(y1,y2),axes=F)
  if ( !any(is.na(shp))) plot( b, add=T, lwd=3)
  rect(x1,x2,y1,(y1+(y2-y1)/3),col="white")
  legend(x=x1,y=(y1+(y2-y1)),fill=rev(col_or),legend=rev(c(br_or[2:length(br_or)])),cex=2.65)
  axis(1,cex.axis=2)
  axis(2,cex.axis=2)
  box()
  text(x=x1+13,y=y2-0.8,labels="(b) 2005-08-28",cex=4)
  par(new=T)
  par(mar=c(3,10,36,3))
  plot(xyFromCell(qqor_a,1:ncell(qqor_a))[,1],getValues(qqor_a),axes=F,pch=21,bg="cornflowerblue",col="cornflowerblue",xlab="",ylab="",
       xlim=c(x1+4.5,x2-2.1),ylim=c(0,300))
  abline(h=seq(0,1000,by=50),col="darkgray",lty=2)
  abline(h=seq(0),col="cornflowerblue",lty=1,lwd=5)
  axis(2,cex.axis=1.5,las=1)
  box()
  dev.off()

  # plot the panel (a) of Fig. 8 in Benestad et al. (2022)
  # precipitation field rescaled according to 1961-1990 climate
  png( file=paste0( "katrina_deflated_1961_1990.png"), height=800, width=1200)
  par(mar=c(3,3,1,1))
  image( qqin1, main="", xlab="", ylab="", col=col_or, breaks=br_or,xlim=c(x1,x2),ylim=c(y1,y2),axes=F)
  if ( !any(is.na(shp))) plot( b, add=T, lwd=3)
  rect(x1,x2,y1,(y1+(y2-y1)/3),col="white")
  box()
  text(x=x1+13,y=y2-0.8,labels="(a) 1961-1990",cex=4)
  par(new=T)
  par(mar=c(3,10,36,3))
  plot(xyFromCell(qqin1_a,1:ncell(qqin1_a))[,1],getValues(qqin1_a),axes=F,pch=21,bg="cornflowerblue",col="cornflowerblue",xlab="",ylab="",
       xlim=c(x1+4.5,x2-2.1),ylim=c(0,300))
  abline(h=seq(0,1000,by=50),col="darkgray",lty=2)
  abline(h=seq(0),col="cornflowerblue",lty=1,lwd=5)
  axis(2,cex.axis=1.5,las=1)
  box()
  dev.off()

  # plot the panel (c) of Fig. 8 in Benestad et al. (2022)
  # precipitation field rescaled to match an hypothetical 2021-2050 climate
  png( file=paste0( "katrina_eninflated_2021_2050.png"), height=800, width=1200)
  par(mar=c(3,3,1,1))
  image( qqin2, main="", xlab="", ylab="", col=col_or, breaks=br_or,xlim=c(x1,x2),ylim=c(y1,y2),axes=F)
  if ( !any(is.na(shp))) plot( b, add=T, lwd=3)
  rect(x1,x2,y1,(y1+(y2-y1)/3),col="white")
  box()
  text(x=x1+13,y=y2-0.8,labels="(c) 2021-2050",cex=4)
  par(new=T)
  par(mar=c(3,10,36,3))
  plot(xyFromCell(qqin2_a,1:ncell(qqin2_a))[,1],getValues(qqin2_a),axes=F,pch=21,bg="cornflowerblue",col="cornflowerblue",xlab="",ylab="",
       xlim=c(x1+4.5,x2-2.1),ylim=c(0,300))
  abline(h=seq(0,1000,by=50),col="gray",lty=2)
  axis(2,cex.axis=1.5,las=1)
  box()
  dev.off()
}

# Function to plot color bar
`color.bar` <- function(col, breaks, nticks=11, title='', cutTails=T,
                      legtxt="",legdig=0,
                      x1=1000000,
                      y1=6450000,
                      x2=1050000,
                      y2=7530000,
                      dx=50000,
                      cex=2.5) {
#    scale = (length(lut)-1)/(max-min)
  nbr<-length(breaks)
  if (cutTails) {
    min<-min(breaks[2:(nbr-1)])
    max<-max(breaks[2:(nbr-1)])
    ticks<-round(seq(2, (nbr-1), len=nticks),0)
  } else {
    min<-min(breaks)
    max<-max(breaks)
    ticks<-round(seq(1, nbr, len=nticks),0)
  }
  dy<-(y2-y1)/length(col)
  rect(x1, y1-dx,
       x2+1.5*dx, y2+dx,
       col="beige", border=NA)
  text((x2+1.5*dx/2),
       (y1+dy/2+(ticks-1)*dy),
       round(breaks[ticks],legdig),
       cex=cex)
  text((x1+x2)/2,
        y2+dx/2,
        legtxt,
        cex=cex)
  for (i in 1:(length(col))) {
   y = (i-1)*dy + y1 
   rect(x1,y,x2,y+dy, col=col[i], border=NA)
  }
#  rect(x1, y1,
#       x2+1.5*dx, y2+dx, border="black")
}

#==============================================================================
# MAIN - MAIN - MAIN - MAIN - MAIN - MAIN - MAIN - MAIN - MAIN - MAIN - MAIN - 
#==============================================================================

#
#------------------------------------------------------------------------------
# Read command line arguments

p <- arg_parser("wave_it")

p <- add_argument(p, "--date1",
                  help="begin (%Y-%m-%d)",
                  type="character",
                  default="2005-08-28")

p <- add_argument(p, "--date2",
                  help="end (%Y-%m-%d)",
                  type="character", 
                  default="2005-08-28")

p <- add_argument(p, "--ffin",
                  help="input file",
                  type="character",
                  default="./data/ERA5_tp_day_20050828.nc")

p <- add_argument(p, "--extent",
                  help="extent",
                  type="numeric",
                  nargs=4,
                  default=c(-107,-70,10,34))

argv <- parse_args(p)

#
#------------------------------------------------------------------------------

# Constants

date_format <- "%Y-%m-%d"
proj4.wgs84 <- "+proj=longlat +datum=WGS84"

#
# time 

if ( !is.na(argv$date1) & !is.na(argv$date2)) {
  Rdate1 <- strptime( paste0(argv$date1, "T00"), format=paste0( date_format, "T%H"), tz="UTC")
  Rdate2 <- strptime( paste0(argv$date2, "T23"), format=paste0( date_format, "T%H"), tz="UTC")
}

#
# open input file

if ( !file.exists( argv$ffin)) q()
nc <-nc_open( argv$ffin)
tp <-nc$var[[2]]
ndims <- tp$ndims
start <- rep(1,ndims)
varsize <- tp$varsize

#
# define timesteps to read from file

tseq <- as.POSIXct( (tp$dim[[3]]$vals * 60 * 60), origin="1900-01-01", tz="UTC")

if ( !is.na(argv$date1) & !is.na(argv$date2)) {
  ixt <- which( as.numeric(tseq) >= as.numeric(Rdate1) & 
                as.numeric(tseq) <= as.numeric(Rdate2))
} else {
  ixt <- 1:length(tseq)
}

#
# define spatial grids

lons_tot <- tp$dim[[1]]$vals
lats_tot <- tp$dim[[2]]$vals
w <- raster( nrows=length(lats_tot), ncols=2*length(lons_tot),  
             xmn=min(lons_tot-360)-0.125, xmx=max(lons_tot)+0.125, 
             ymn=min(lats_tot)-0.125, ymx=max(lats_tot)+0.125, 
             crs=proj4.wgs84)
r <- crop(w, extent(c(-180.1, 179.8, -90.125, 90.125)))
xy<-xyFromCell(r,1:ncell(r))
lons<-sort(unique(xy[,1]))
lats<-sort(unique(xy[,2]),decreasing=T)
rm(xy)

# number of rows (cols) of the dyadic grid 
dimdy  <- 2**ceiling( log2( max( dim(r)[1:2])))

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
# elaboration

nt      <- length(ixt)
En2o_t<-array(data=NA,dim=c(nnscales,nt))

# loop over the time steps 
for( i in 1:nt ) {
  if ( ( i %% 10) == 0) cat( paste(i,"/",nt,"\n"))
  # read gridded data from nc-file
  start <- rep(1,ndims)
  start[ndims] <- ixt[i]
  count <- varsize
  count[ndims] <- 1
  data <- ncvar_get( nc, tp, start=start, count=count )
  # trick to get the precipitation field onto the expected grid
  w[] <- cbind( t(data), t(data)) * 1000
  r <- crop(w, extent(c(-180.1, 179.8, -90.125, 90.125)))
  # resample the precipitation field onto the dyadic domain
  u <- resample( r, s, method="bilinear")
  # store the raster on a dedicated environment
  env <- new.env( parent = emptyenv())
  env$r <- u
  # call plotting function
#  plot_multires(env,shp=c("/home/cristianl/data/geoinfo/TM_WORLD_BORDERS_LATLON/TM_WORLD_BORDERS-0.2.shp","TM_WORLD_BORDERS-0.2"))
  plot_multires(env)
} # end loop over the time steps

#
# close the input file
 
nc_close(nc)

#
# Normal exit

q( status=0)
