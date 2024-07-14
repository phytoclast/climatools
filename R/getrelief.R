# library(climatools)
# library(terra)
# dem <- rast('C:/a/geo/dem/SRTM250m/usa/w001001.adf')
#
# dem1 <-  reproject(dem=dem, rs = 250, w=500000, h=500000, lat=47, lon=-122)
# plot(dem1)
#
# dem2 <-  reproject(dem=dem, rs = 250, w=500000, h=500000, lat=44, lon=-85)
# plot(dem2)


#' Get topographic relief
#'
#' @param dm Digital elevation model (terra raster) with projected coordinates.
#' @param r1 Minimum focal radius for estimating relief (only focal radius when n-0).
#' @param r2 Maximum focal radius for estimating relief.
#' @param n Number of segments between maximum and minimum radii to estimate relief.
#' @param s Slope for which relief is estimated (Default is 0.1 for a 10% slope). Should be a ratio of vertical units to horizontal units).
#'
#' @return Topographic relief of either a fixed radius or for a fixed slope for a range of radii.
#' @export
#'
#' @examples
getrelief <- function(dm, r1, r2, n=0, s=0.1){
  require(terra)
  if(n==0){#simple relief using single fixed radius
    xmax <- focalmax(dm, r1/2)
    xmin <- focalmin(dm, r1/2)
    rng <- xmax - xmin
  }else{#relief maintaining a slope within a range of radii
    d = log(r2) -  log(r1)
    ser <- n:0
    rx <- exp(log(r1)+d/n*ser)
    for(i in 1:length(rx)){#i=2
      xmax <- focalmax(dm, rx[i]/2)
      xmin <- focalmin(dm, rx[i]/2)
      if(i==1){
        rng1 <- xmax - xmin
        rx1 <- rx[i]
        rate1 <- (rng1/(rx1*s))
        rng0 <- rng1*(rate1>=1)
      }else{
        rng2 <- rng1
        rng1 <- xmax - xmin
        rx2 <- rx1
        rx1 <- rx[i]
        rate2 <- (rng2/(rx2*s))
        rate1 <- (rng1/(rx1*s))

        num <- (rate1-1)
        den <- (rate1-rate2)

        den <- ifel(den >0,den,0)
        numden <- num/(den+0.00001)
        if(i!=length(rx)){
          rng <- (numden*rng2 + (1-numden)*rng1)*(rate1>=1)*(rng0<=0)+rng0*(rng0>0)
          rng0 <- rng
        }else{
          rng <- (rate1<1)*(rng0<=0)*rng1+(numden*rng2 + (1-numden)*rng1)*(rate1>=1)*(rng0<=0)+rng0*(rng0>0)}
      }}}
  return(rng)}
