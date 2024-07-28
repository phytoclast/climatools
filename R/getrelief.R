#' Get topographic relief
#'
#' @param dm Digital elevation model (terra raster) with projected coordinates.
#' @param r1 Minimum focal radius for estimating relief (only focal radius when n-0).
#' @param r2 Maximum focal radius for estimating relief.
#' @param n Number of segments between maximum and minimum radii to estimate relief.
#' @param s Slope for which relief is estimated (Default is 0.1 for a 10\% slope). Should be a ratio of vertical units to horizontal units).
#' @param p Precision from low to exact with higher levels of precision requiring more processing time.
#' @param breaks Optional vector of relief thresholds to specifically calculate a radius rather than to depend on interpolation among a regular interval between maximum and minimum radii.
#'
#' @return Topographic relief of either a fixed radius or for a fixed slope for a range of radii.
#'
#' @export
#'
#' @examples library(terra)
#' #create fake raster
#' x <- (1:10)/10-0.5
#' y <- x
#' df <- merge(x,y)
#' df$z <- 1000/((df$x-0.25)^2+(df$y-0.40)^2+1)^2
#' set.seed(1)
#' df$z <- df$z +rnorm(10*10, sd = 100)
#' df$x <- df$x*10 -85
#' df$y <- df$y*10 +20
#' dem <- rast(cbind(x=df$x,y=df$y,z=df$z), type="xyz", crs='EPSG:4326')
#' dem <- reproject(dem=dem, rs = 250)
#' plot(dem)
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
getrelief <- function(dm, r1, r2, n=0, s=0.1, p=c('low', 'medium', 'high','exact'), breaks=NA){
  require(terra)
  #p is for precision options
  p=p[1]

  if(n==0){#simple relief using single fixed radius
    xmax <- focalmax(dm, r1/2,p=p)
    xmin <- focalmin(dm, r1/2,p=p)
    rng <- xmax - xmin
  }else{#relief maintaining a slope within a range of radii
    d = log(r2) -  log(r1)
    ser <- n:0
    rx <- exp(log(r1)+d/n*ser)
    if(!is.na(breaks[1])){
    rx <- unique(sort(c(breaks/s,rx), decreasing=T))}

    for(i in 1:length(rx)){#i=2
      xmax <- focalmax(dm, rx[i]/2, p=p)
      xmin <- focalmin(dm, rx[i]/2, p=p)
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
