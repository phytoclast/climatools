#' Get topographic relief
#'
#' @param dm Digital elevation model (terra raster) with projected coordinates.
#' @param r1 Minimum focal radius for estimating relief (only focal radius when n-0). Optional if using breaks.
#' @param r2 Maximum focal radius for estimating relief. Optional if using breaks.
#' @param n Number of segments between maximum and minimum radii to estimate relief. n=0 for a simple focal range analysis using r1.
#' @param s Slope for which relief is estimated (Default is 0.1 for a 10\% slope). Should be a ratio of vertical units to horizontal units).
#' @param p Precision from low to exact with higher levels of precision requiring more processing time.
#' @param breaks Optional vector of relief thresholds to specifically calculate a radius rather than to depend on interpolation among a regular interval between maximum and minimum radii. Will combine with r1 and r2 vertical equivalence if they are populated.
#' @param reverse If true, give priority to widest focal neighborhoods, otherwise relief has to meet slope criteria for all the smaller focal neigborhoods first.
#'
#' @return Topographic relief of either a fixed radius or for a fixed slope for a range of radii. Relief above and below focal radius breaks is simple, not based on fixed slope.
#'
#' @export
#'
#' @examples library(terra)#'
#' #load package data for elevation
#' data("denali")
#' dem <- denali
#' dem <- toraster(dem)
#' plot(dem)
#' #Reproject to get metric distance units.
#' dem <- reproject(dem, rs=250)
#'
#' #Established desired vertical increment of focal range analysis.
#' prebreaks = c(50, 100,150,200,300,500,1000,2000,3000,4000,5000,6000,7000,8000)
#' maxrange <- max(minmax(dem)[2],1000)#Establish maximum relative to the highest elevation in DEM. This is to reduce the need for radii wider than required.
#' breaks <- prebreaks[prebreaks <=maxrange]
#' #Relief at a 10% slope. Using "low" precision for faster performance
#' rng10 <- getrelief(dem, r1=1000, r2=maxrange*10, s=0.1, n=1, p='low', breaks = breaks)
#' plot(rng10)
#' #Relief at a 25% slope. Using "medium" precision for better quality.
#' rng25 <- getrelief(dem, r1=1000, r2=maxrange*10, s=0.25, n=1, p='medium', breaks = breaks)
#' plot(rng25)
#'
#' #get relief with a 10% slope using 10 intervals between focal radii and default low precision.
#' rng <- getrelief(dem, r1=1000, r2=100000, s=0.1, n=10)
#' #get relief with a 25% slope using 10 intervals between focal radii and default low precision.
#' rng <- getrelief(dem, r1=1000, r2=100000, s=0.25, n=10)
#' plot(rng)
#' #get 10% relief using specific break points of vertical relief and medium precision.
#' rng <- getrelief(dem, r1=1000, r2=100000, s=0.1, n=1, p='medium', breaks = c(300, 1000, 2500, 5000))
#' plot(rng)
#' #get relief with a fixed radius of 2000 meters.
#' rng <- getrelief(dem, r1=2000, s=0.1, n=0)
#' plot(rng)
#'
getrelief <- function(dm, r1=NA, r2=NA, n=0, s=0.1, p=c('low', 'medium', 'high','exact'), breaks=NA, reverse = FALSE){
  require(terra)
  #p is for precision options
  p=p[1]

  if(n==0){#simple relief using single fixed radius
    xmax <- focalmax(dm, r1/2,p=p)
    xmin <- focalmin(dm, r1/2,p=p)
    rng <- xmax - xmin
  }else{#relief maintaining a slope within a range of radii
    if(!is.na(r1)&!is.na(r2)&!is.na(breaks[1])){
      d = log(r2) -  log(r1)
      ser <- n:0
      rx <- exp(log(r1)+d/n*ser)
      rx <- unique(round(sort(c(breaks/s,rx), decreasing=reverse),1))
    }else if(!is.na(breaks[1])){
      #determine maximum range of input raster to save processing time
      rnmax <- terra::minmax(dm)[2]-terra::minmax(dm)[1]
      rx <- unique(round(sort(c(breaks/s,rnmax/s), decreasing=reverse),1)); rx <- rx[rx <= rnmax/s]
    }else if(!is.na(r1)&!is.na(r2)){d = log(r2) -  log(r1)
    ser <- n:0
    rx <- exp(log(r1)+d/n*ser)}else{message("Not enough parameters.")}
    if(reverse){
      for(i in 1:length(rx)){#i=1
        xmax <- focalmax(dm, rx[i]/2, p=p)
        xmin <- focalmin(dm, rx[i]/2, p=p)
        if(i==1){#i=1
          rng <- xmax - xmin
          rx0 <- rx[i]*s
        }else{#i=3
          rng2 <- rng
          rng1 <- xmax - xmin
          rx2 <- rx0
          rx1 <- rx[i]*s
          rx0 <- rx1
          wt0 <- (rng1 < rx1)*(rng2 < rx2)
          wt1 <- (rx2-rng2)*(rng2 < rx2)*(rng1 >= rx1)
          wt2 <- (rng1-rx1)*(rng2 < rx2)*(rng1 >= rx1)
          wt3 <- (rng2 >= rx2)
          rng <- (rng1*wt0+rng1*wt1+rng2*wt2+rng2*wt3)/(wt0+wt1+wt2+wt3)
        }}
    }else{
      for(i in 1:length(rx)){#i=1
        xmax <- focalmax(dm, rx[i]/2, p=p)
        xmin <- focalmin(dm, rx[i]/2, p=p)
        if(i==1){#i=1
          rng <- xmax - xmin
          rx0 <- rx[i]*s
        }else{#i=3
          rng2 <- xmax - xmin
          rng1 <- rng
          rx2 <- rx[i]*s
          rx1 <- rx0
          rx0 <- rx2
          wt0 <- (rng1 < rx1)
          wt1 <- (rx2-rng2)*(rng2 < rx2)*(rng1 >= rx1)
          wt2 <- (rng1-rx1)*(rng2 < rx2)*(rng1 >= rx1)
          wt3 <- (rng2 >= rx2)*(rng1 >= rx1)
          rng <- (rng1*wt0+rng1*wt1+rng2*wt2+rng2*wt3)/(wt0+wt1+wt2+wt3)
        }}
    }

  }
  return(rng)}

# original_getrelief <- function(dm, r1, r2, n=0, s=0.1, p=c('low', 'medium', 'high','exact'), breaks=NA){
#   require(terra)
#   #p is for precision options
#   p=p[1]
#
#   if(n==0){#simple relief using single fixed radius
#     xmax <- focalmax(dm, r1/2,p=p)
#     xmin <- focalmin(dm, r1/2,p=p)
#     rng <- xmax - xmin
#   }else{#relief maintaining a slope within a range of radii
#     d = log(r2) -  log(r1)
#     ser <- n:0
#     rx <- exp(log(r1)+d/n*ser)
#     if(!is.na(breaks[1])){
#       rx <- unique(sort(c(breaks/s,rx), decreasing=T))}
#
#     for(i in 1:length(rx)){#i=2
#       xmax <- focalmax(dm, rx[i]/2, p=p)
#       xmin <- focalmin(dm, rx[i]/2, p=p)
#       if(i==1){
#         rng1 <- xmax - xmin
#         rx1 <- rx[i]
#         rate1 <- (rng1/(rx1*s))
#         rng0 <- rng1*(rate1>=1)
#       }else{
#         rng2 <- rng1
#         rng1 <- xmax - xmin
#         rx2 <- rx1
#         rx1 <- rx[i]
#         rate2 <- (rng2/(rx2*s))
#         rate1 <- (rng1/(rx1*s))
#
#         num <- (rate1-1)
#         den <- (rate1-rate2)
#
#         den <- ifel(den >0,den,0)
#         numden <- num/(den+0.00001)
#         if(i!=length(rx)){
#           rng <- (numden*rng2 + (1-numden)*rng1)*(rate1>=1)*(rng0<=0)+rng0*(rng0>0)
#           rng0 <- rng
#         }else{
#           rng <- (rate1<1)*(rng0<=0)*rng1+(numden*rng2 + (1-numden)*rng1)*(rate1>=1)*(rng0<=0)+rng0*(rng0>0)}
#       }}}
#   return(rng)}
#
# ####This version allows neigborhood to be redistributed to emphasize either upslope areas versus downslope areas.
# x_getrelief <- function(dm, r1=NA, r2=NA, n=0, s=0.1, p=c('low', 'medium', 'high','exact'), breaks=NA, reverse = FALSE, upratio = 0.5){
#   require(terra)
#   #p is for precision options
#   p=p[1]
#   if(upratio > 1 | upratio <0){message("Upratio needs to be between 0 and 1")}
#
#   if(n==0){#simple relief using single fixed radius
#     xmax <- focalmax(dm, r1*upratio,p=p)
#     xmin <- focalmin(dm, r1*(1-upratio),p=p)
#     rng <- xmax - xmin
#   }else{#relief maintaining a slope within a range of radii
#     if(!is.na(r1)&!is.na(r2)&!is.na(breaks[1])){
#       d = log(r2) -  log(r1)
#       ser <- n:0
#       rx <- exp(log(r1)+d/n*ser)
#       rx <- unique(round(sort(c(breaks/s,rx), decreasing=reverse),1))
#     }else if(!is.na(breaks[1])){
#       #determine maximum range of input raster to save processing time
#       rnmax <- terra::minmax(dm)[2]-terra::minmax(dm)[1]
#       rx <- unique(round(sort(c(breaks/s,rnmax/s), decreasing=reverse),1)); rx <- rx[rx <= rnmax/s]
#     }else if(!is.na(r1)&!is.na(r2)){d = log(r2) -  log(r1)
#     ser <- n:0
#     rx <- exp(log(r1)+d/n*ser)}else{message("Not enough parameters.")}
#     if(reverse){
#       for(i in 1:length(rx)){#i=1
#         xmax <- focalmax(dm, rx[i]*upratio, p=p)
#         xmin <- focalmin(dm, rx[i]*(1-upratio)/2, p=p)
#         if(i==1){#i=1
#           rng <- xmax - xmin
#           rx0 <- rx[i]*s
#         }else{#i=3
#           rng2 <- rng
#           rng1 <- xmax - xmin
#           rx2 <- rx0
#           rx1 <- rx[i]*s
#           rx0 <- rx1
#           wt0 <- (rng1 < rx1)*(rng2 < rx2)
#           wt1 <- (rx2-rng2)*(rng2 < rx2)*(rng1 >= rx1)
#           wt2 <- (rng1-rx1)*(rng2 < rx2)*(rng1 >= rx1)
#           wt3 <- (rng2 >= rx2)
#           rng <- (rng1*wt0+rng1*wt1+rng2*wt2+rng2*wt3)/(wt0+wt1+wt2+wt3)
#         }}
#     }else{
#       for(i in 1:length(rx)){#i=1
#         xmax <- focalmax(dm, rx[i]*upratio, p=p)
#         xmin <- focalmin(dm, rx[i]*(1-upratio)/2, p=p)
#         if(i==1){#i=1
#           rng <- xmax - xmin
#           rx0 <- rx[i]*s
#         }else{#i=3
#           rng2 <- xmax - xmin
#           rng1 <- rng
#           rx2 <- rx[i]*s
#           rx1 <- rx0
#           rx0 <- rx2
#           wt0 <- (rng1 < rx1)
#           wt1 <- (rx2-rng2)*(rng2 < rx2)*(rng1 >= rx1)
#           wt2 <- (rng1-rx1)*(rng2 < rx2)*(rng1 >= rx1)
#           wt3 <- (rng2 >= rx2)*(rng1 >= rx1)
#           rng <- (rng1*wt0+rng1*wt1+rng2*wt2+rng2*wt3)/(wt0+wt1+wt2+wt3)
#         }}
#     }
#
#   }
#   return(rng)}
