#Digital elevation model focal statistics tools using terra package
#


#focal neighborhood circle smoother than terra::focalMat() function
focalCircle <- function(x, r){
  rs <- res(x)[1]
  n = floor(r/rs)*2+1
  mx <- matrix(0, nrow = n, ncol = n)
  for(i in 1:n){
    for(k in 1:n){
      k0 <- (k - ceiling(n/2))/n*2
      i0 <- (i - ceiling(n/2))/n*2
      mx[i,k] <- (k0^2+i0^2) <= 1
    }
  }
  return(mx)}










#' Get focal maximum with specified radius
#'
#' @param x raster
#' @param r radius
#' @param p precision, from low to exact, with higher levels of precision requiring more processing time.
#'
#' @return focal maximum raster
#' This function makes use of lower resolution aggregated rasters to speed up calculations at increasing radii, compromising accuracy. Higher precision generates a focal neighborhood at a higher resolution, while lower precision aggregates more aggressively before focal analysis.
#' @export
#'
#' @examples
focalmax <- function(x, r, p=c('low', 'medium', 'high','exact')){
  require(terra)
  #establish aggregating factor when radius is too large
  #p is for precision options
  p=p[1]
  if(p == 'low'){
    fc = floor(r/res(x)[1]/9+1)
  }else if(p == 'medium'){
    fc = floor(r/res(x)[1]/21+1)
  }else if(p == 'high'){
    fc = floor(r/res(x)[1]/81+1)
  }else{
    fc = floor(r/res(x)[1]/243+1)}

  x1 = x
  if(fc > 1){
    x1 <- aggregate(x, fact = fc, fun = 'max',  na.rm=TRUE)
  }
  #create a focal weights matrix of the appropriate size for given radius
  #fm <- focalMat(x1, d=r, type = 'circle')
  fm <- focalCircle(x1, r=r)
    #exclude outer portion of circle and ensure max/min values are only multplied by 1
  fm.na <- ifelse(fm > 0, 1, NA)
  #reduce mat size if raster too small
  matsize = nrow(fm.na)
  centermat <- floor(nrow(fm.na)/2+1)
  goodrows <- intersect((1:matsize),(centermat-floor(nrow(x1)-1)):(centermat+floor(nrow(x1)-1)))
  goodcols <- intersect((1:matsize),(centermat-floor(ncol(x1)-1)):(centermat+floor(ncol(x1)-1)))
  fm.na <- fm.na[goodrows,goodcols, drop = FALSE]
  if(nrow(fm.na) <= 1 & ncol(fm.na) <= 1){
    x1.max <- x1
  }else{
    x1.max <- focal(x1, fm.na, fun='max', na.rm=T)}
  #restore resolution in result
  if(fc > 1){
    x1.max <- project(x1.max, x)
  }
  return(x1.max)}

#' Get focal minimum with specified radius
#'
#' @param x raster
#' @param r radius
#' @param p precision, from low to exact, with higher levels of precision requiring more processing time.
#'
#' @return focal minimum raster
#' This function makes use of lower resolution aggregated rasters to speed up calculations at increasing radii, compromising accuracy. Higher precision generates a focal neighborhood at a higher resolution, while lower precision aggregates more aggressively before focal analysis.
#' @export
#'
#' @examples
focalmin <- function(x, r, p=c('low', 'medium', 'high','exact')){
  require(terra)
  #establish aggregating factor when radius is too large
  #p is for precision options
  p=p[1]
  if(p == 'low'){
    fc = floor(r/res(x)[1]/9+1)
  }else if(p == 'medium'){
    fc = floor(r/res(x)[1]/21+1)
  }else if(p == 'high'){
    fc = floor(r/res(x)[1]/81+1)
  }else{
    fc = floor(r/res(x)[1]/243+1)}

  x1 = x
  if(fc > 1){
    x1 <- aggregate(x, fact = fc, fun = 'min',  na.rm=TRUE)
  }
  #create a focal weights matrix of the appropriate size for given radius
  #fm <- focalMat(x1, d=r, type = 'circle')
  fm <- focalCircle(x1, r=r)
  #exclude outer portion of circle and ensure max/min values are only multplied by 1
  fm.na <- ifelse(fm > 0, 1, NA)
  #reduce mat size if raster too small
  matsize = nrow(fm.na)
  centermat <- floor(nrow(fm.na)/2+1)
  goodrows <- intersect((1:matsize),(centermat-floor(nrow(x1)-1)):(centermat+floor(nrow(x1)-1)))
  goodcols <- intersect((1:matsize),(centermat-floor(ncol(x1)-1)):(centermat+floor(ncol(x1)-1)))
  fm.na <- fm.na[goodrows,goodcols, drop = FALSE]
  if(nrow(fm.na) <= 1 & ncol(fm.na) <= 1){
    x1.min <- x1
  }else{
    x1.min <- focal(x1, fm.na, fun='min', na.rm=T)}
  #restore resolution in result
  if(fc > 1){
    x1.min <- project(x1.min, x)
  }
  return(x1.min)}



#' Get focal median with specified radius
#'
#' @param x raster
#' @param r radius
#' @param p precision, from low to exact, with higher levels of precision requiring more processing time.
#'
#' @return focal median raster
#' This function makes use of lower resolution aggregated rasters to speed up calculations at increasing radii, compromising accuracy. Higher precision generates a focal neighborhood at a higher resolution, while lower precision aggregates more aggressively before focal analysis.
#' @export
#'
#' @examples
focalmed <- function(x, r, p=c('low', 'medium', 'high','exact')){
  require(terra)
  #establish aggregating factor when radius is too large
  #p is for precision options
  p=p[1]
  if(p == 'low'){
    fc = floor(r/res(x)[1]/9+1)
  }else if(p == 'medium'){
    fc = floor(r/res(x)[1]/21+1)
  }else if(p == 'high'){
    fc = floor(r/res(x)[1]/81+1)
  }else{
    fc = floor(r/res(x)[1]/243+1)}

  x1 = x
  if(fc > 1){
    x1 <- aggregate(x, fact = fc, fun = 'mean',  na.rm=TRUE)
  }
  #create a focal weights matrix of the appropriate size for given radius
  #fm <- focalMat(x1, d=r, type = 'circle')
  fm <- focalCircle(x1, r=r)
  #exclude outer portion of circle and ensure max/min values are only multplied by 1
  fm.na <- ifelse(fm > 0, 1, NA)
  #reduce mat size if raster too small
  matsize = nrow(fm.na)
  centermat <- floor(nrow(fm.na)/2+1)
  goodrows <- intersect((1:matsize),(centermat-floor(nrow(x1)-1)):(centermat+floor(nrow(x1)-1)))
  goodcols <- intersect((1:matsize),(centermat-floor(ncol(x1)-1)):(centermat+floor(ncol(x1)-1)))
  fm.na <- fm.na[goodrows,goodcols, drop = FALSE]
  if(nrow(fm.na) <= 1 & ncol(fm.na) <= 1){
    x1.mean <- x1
  }else{
    x1.mean <- focal(x1, fm.na, fun='mean', na.rm=T)}
  #restore resolution in result
  if(fc > 1){
    x1.mean <- project(x1.mean, x)
  }
  return(x1.mean)}


#' Calculate hillslope position
#'
#' @param dm original raster
#' @param r neighborhood radius
#'
#' @return Raster with relative slope position 0 to 1.
#' @export
#'
#' @examples
hillpos <- function(dm, r){#relative slope position
  xmax = focalmax(dm, r)
  xmin = focalmin(dm, r)
  xmed = focalmed(dm, r)
  x.pos <- (dm - xmed)/(xmax - xmed+0.5)
  x.neg <- (dm - xmed)/(xmed - xmin+0.5)
  x.pos <- ifel(x.pos > 0, x.pos,0)
  x.neg <- ifel(x.neg < 0, x.neg,0)
  p <- ((x.pos+x.neg)+1)/2
  return(p)
}

#compound slope position
#' Compound hillslope position using 3 neighborhood scales
#'
#' @param dm original raster
#' @param r1 small neighborhood radius
#' @param r2 intermediate neighborhood radius
#' @param r3 large neighborhood radius
#'
#' @return Compound hillslope position raster
#' This function gives more weight to neighborhoods with higher relative relief.
#' @export
#'
#' @examples
comphillpos = function(dm, r1, r2, r3){

  x.pos1 <- hillpos(dm, r1)
  x.pos2 <- hillpos(dm, r2)
  x.pos3 <- hillpos(dm, r3)
  x.pos.r1 <- focalmax(x.pos2, r1) - focalmin(x.pos2, r1)
  x.pos.1 <- x.pos1*x.pos.r1 + x.pos2*(x.pos.r1*-1+1)
  x.pos.r2 <- focalmax(x.pos3, r2) - focalmin(x.pos3, r2)
  x.pos <- x.pos.1*x.pos.r2 + x.pos3*(x.pos.r2*-1+1)
  return(x.pos)
}

#' Topographic position index using 3 neighborhood scales
#'
#' @param dm original raster
#' @param r1 small neighborhood radius
#' @param r2 intermediate neighborhood radius
#' @param r3 large neighborhood radius
#'
#' @return Topographic position index raster
#' This function simply averages the hillslope positions for all 3 neighborhoods.
#' @export
#'
#' @examples
tpi = function(dm, r1, r2, r3){
  x.pos1 <- hillpos(dm, r1)
  x.pos2 <- hillpos(dm, r2)
  x.pos3 <- hillpos(dm, r3)
  x.pos <- (x.pos1+x.pos2+x.pos3)/3
  return(x.pos)
}

#enhance raster is designed to increase precision
#t1 = climate grid (mainly temperature as precipitation has a fuzzier relationship with elevation)

#e1 = elevation grid matching resolution of climate grid

#e2 = elevation grid of higher resolution and cropped to an area of interest

#' Enhance raster using trends with a higher resolution digital elevation model.
#'
#' @param t1 Climate or other types of grids that should have a tight relationship with elevation.
#' @param e1 Elevation grid matching resolution and extent of climate grid of interest.
#' @param e2 Elevation grid of higher resolution, cropped to an area of interest.
#'
#' @return High resolution climate grid.
#' @export
#'
#' @examples
enhanceRast <- function(t1,e1,e2){
  #create new extent to crop analysis
  expts <- data.frame(x = c(ext(e2)[1],ext(e2)[1],ext(e2)[2],ext(e2)[2]), y = c(ext(e2)[3],ext(e2)[4],ext(e2)[3],ext(e2)[4]))
  #extend a little
  expts <- expts + matrix(c(-10000,-10000,10000,10000,-10000,10000,-10000,10000),ncol = 2)
  #convert to spatVect and project to get new extent
  expts <- vect(expts, geom=c('x','y'), crs=crs(e2))
  expts <- project(expts, e1)
  #apply new extent to crop analysis
  e1 <- crop(e1, ext(expts))
  t1 <- crop(t1, ext(expts))
  #filter out bogus elevations
  e1[e1 > 9000] <- NA; e1[e1 < -500] <- NA
  #ensure grids match extent and resolution
  e1 <- project(e1, t1)
  #ensure that grid is numeric not factor
  e1 <- e1+0

  names(t1) <- 't1';names(e1) <- 'e1';names(e2) <- 'e2'
  #find neighborhood means
  t5km <- aggregate(t1, fact=10, fun='mean', na.rm=T)
  e5km <- aggregate(e1, fact=10, fun='mean', na.rm=T)
  t5km <- resample(t5km, t1, method='near'); names(t5km)<- 't5km'
  e5km <- resample(e5km, t1, method='near'); names(e5km)<- 'e5km'
  #difference from neighborhood means
  tdif <- t1-t5km
  edif <- e1-e5km
  #average lapse rates per elevation in a aggregate neighborhood with higher weights where elevation differences are greatest, then smoothed and resampled to target resolution
  wts <- edif^2
  rate <- tdif/(edif+0.1)
  rate.wts <- wts*rate
  rate.sum <- aggregate(rate.wts, fact=10, fun='sum', na.rm=T)
  wts.sum <- aggregate(wts, fact=10, fun='sum', na.rm=T)
  rate.5km <- rate.sum/(wts.sum+0.001)
  rate.5km <- focal(rate.5km, na.rm=T, fun="median")
  t1.90 <- project(t1, e2)
  e1.90 <- project(e1, e2)
  rate.90 <- project(rate.5km, e2)
  #apply lapse rate to high resolution DEM
  new.90 <- t1.90 + (e2 - e1.90)*rate.90
  return(new.90)}


#Reduces the resolution of a raster but preserves the values of local high and low values.
#' Reduction in horizontal resolution preserving vertical range.
#'
#' @param hires Input high resolution raster.
#' @param fact Aggregation factor to reduce resolution by.
#'
#' @return Reduced resolution raster maintaining the vertical range of the parent raster.
#' @export
#'
#' @examples
AmplifiedReduction <- function(hires, fact = 5){

  lowres <- aggregate(hires, fact=fact, fun='mean',na.rm=TRUE)

  hiresMax <- aggregate(hires, fact=fact, fun='max',na.rm=TRUE)
  hiresMin <- aggregate(hires, fact=fact, fun='min',na.rm=TRUE)

  lowresMax <- aggregate(lowres, fact=fact, fun='max',na.rm=TRUE)
  lowresMin <- aggregate(lowres, fact=fact, fun='min',na.rm=TRUE)
  lowresMean <- aggregate(lowres, fact=fact, fun='mean',na.rm=TRUE)

  lowresMax <- resample(lowresMax, lowres, method = 'bilinear')
  lowresMin <- resample(lowresMin, lowres, method = 'bilinear')
  lowresMean <- resample(lowresMean, lowres, method = 'bilinear')

  ElevBin <- lowres >= lowresMean

  ElevRel1 <- (lowres - lowresMean)/(lowresMax - lowresMin + 0.1)
  ElevRel2 <- (lowresMean - lowres)/(lowresMean - lowresMin + 0.1)

  ElevAmp1 <- (hiresMax * (ElevRel1) + lowres * (ElevRel1-1)*-1)
  ElevAmp2 <- (lowres * (ElevRel2) + hiresMin * (ElevRel2-1)*-1)

  ElevAmp <- ElevAmp1 * ElevBin + ElevAmp2 * (ElevBin-1)*-1
  return(ElevAmp)}

# EnhancedResample <- function(e1,e2){
#   #probably create maximum and minimum rasters of original raster, and then a relative elevation raster from the resampled/reprojected raster.
# }
