#Digital elevation model focal statistics tools using terra package
#


#focal neighborhood circle smoother than terra::focalMat() function
focalCircle <- function(x, r){
  u <- linearUnits(x)
  u <- ifelse(u == 0, 10000000/90,u)
  rs <- res(x)[1]*u
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
#' @param r radius (meters)
#' @param p precision, from low to exact, with higher levels of precision requiring more processing time.
#'
#' @return focal maximum raster
#' This function makes use of lower resolution aggregated rasters to speed up calculations at increasing radii, compromising accuracy. Higher precision generates a focal neighborhood at a higher resolution, while lower precision aggregates more aggressively before focal analysis.
#' @export
#'
#' @examples
#' data("denali")
#' dem <- denali
#' dem <- toraster(dem)
#' dem <- reproject(dem, rs=250)
#' plot(dem)
#' x <- climatools::focalmax(dem, 2000)
#' plot(x)
focalmax <- function(x, r, p=c('low', 'medium', 'high','exact')){
  require(terra)
  #establish aggregating factor when radius is too large
  #p is for precision options
  u <- linearUnits(x)
  u <- ifelse(u == 0, 10000000/90,u)
  rs <- res(x)[1]*u
  p=p[1]
  if(p == 'low'){
    fc = floor(r/rs/9+1)
  }else if(p == 'medium'){
    fc = floor(r/rs/21+1)
  }else if(p == 'high'){
    fc = floor(r/rs/81+1)
  }else{
    fc = floor(r/rs/243+1)}

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
#' @param r radius (meters)
#' @param p precision, from low to exact, with higher levels of precision requiring more processing time.
#'
#' @return focal minimum raster
#' This function makes use of lower resolution aggregated rasters to speed up calculations at increasing radii, compromising accuracy. Higher precision generates a focal neighborhood at a higher resolution, while lower precision aggregates more aggressively before focal analysis.
#' @export
#'
#' @examples
#' data("denali")
#' dem <- denali
#' dem <- toraster(dem)
#' dem <- reproject(dem, rs=250)
#' plot(dem)
#' x <- climatools::focalmin(dem, 2000)
#' plot(x)
focalmin <- function(x, r, p=c('low', 'medium', 'high','exact')){
  require(terra)
  #establish aggregating factor when radius is too large
  #p is for precision options
  u <- linearUnits(x)
  u <- ifelse(u == 0, 10000000/90,u)
  rs <- res(x)[1]*u
  p=p[1]
  if(p == 'low'){
    fc = floor(r/rs/9+1)
  }else if(p == 'medium'){
    fc = floor(r/rs/21+1)
  }else if(p == 'high'){
    fc = floor(r/rs/81+1)
  }else{
    fc = floor(r/rs/243+1)}

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
#' @param r radius (meters)
#' @param p precision, from low to exact, with higher levels of precision requiring more processing time.
#'
#' @return focal median raster
#' This function makes use of lower resolution aggregated rasters to speed up calculations at increasing radii, compromising accuracy. Higher precision generates a focal neighborhood at a higher resolution, while lower precision aggregates more aggressively before focal analysis.
#' @export
#'
#' @examples
#' data("denali")
#' dem <- denali
#' dem <- toraster(dem)
#' dem <- reproject(dem, rs=250)
#' plot(dem)
#' x <- climatools::focalmed(dem, 2000)
#' plot(x)
focalmed <- function(x, r, p=c('low', 'medium', 'high','exact')){
  require(terra)
  #establish aggregating factor when radius is too large
  #p is for precision options
  u <- linearUnits(x)
  u <- ifelse(u == 0, 10000000/90,u)
  rs <- res(x)[1]*u
  p=p[1]
  if(p == 'low'){
    fc = floor(r/rs/9+1)
  }else if(p == 'medium'){
    fc = floor(r/rs/21+1)
  }else if(p == 'high'){
    fc = floor(r/rs/81+1)
  }else{
    fc = floor(r/rs/243+1)}

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
#' data("denali")
#' dem <- denali
#' dem <- toraster(dem)
#' dem <- reproject(dem, rs=250)
#' plot(dem)
#' x <- climatools::hillpos(dem, 2000)
#' plot(x)
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
#' data("denali")
#' dem <- denali
#' dem <- toraster(dem)
#' dem <- reproject(dem, rs=250)
#' plot(dem)
#' x <- climatools::comphillpos(dem, 500, 2000, 5000)
#' plot(x)
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
#' data("denali")
#' dem <- denali
#' dem <- toraster(dem)
#' dem <- reproject(dem, rs=250)
#' plot(dem)
#' x <- climatools::tpi(dem, 500, 2000, 5000)
#' plot(x)
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
#' @examples #load package data for elevation
#' data("denali")
#' dem0 <- denali
#' dem0 <- toraster(dem0)
#' dem <- reproject(dem0, rs=250)
#' plot(dem)
#'
#' #load package data for July temperature
#' data("Tw")
#' temperature0 <- Tw
#' temperature0 <- toraster(temperature0)
#'
#' #Original temperature is coarse resolution.
#' plot(temperature0)
#'
#' #This resamples to match resolution of elevation, but it appears blurry.
#' temperature <- project(temperature0, dem)
#' plot(temperature)
#'
#' #prepare lower resolution elevation raster matching temperature resolution
#' dem0 <- project(dem0, temperature0)
#' plot(dem0)

#' #This estimates a local relationship of temperature with elevation.
#' temperature <- enhanceRast(temperature0, dem0,  dem)
#' plot(temperature)
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

#' Restore maximum and minimum values of a resampled raster
#'
#' @param x Resampled or reprojected raster (lost resolution).
#' @param y Original raster (high resolution).
#' @param s Optional point simple feature vector (sf) with locations of known high points, to further enhance accuracy of elevation extremes that may have been smoothed out in DEM resampling.
#' @param e Optional name of the column corresponding to elevation in the simple feature vector.
#'
#' @return Resampled raster with original maximum and minimum extreme values restored.
#' This function is used to correct for the loss in highest and lowest values within a local neighborhood after resampling (i.e. bilinear method) blends them with adjacent values (or substitutes with nearest neighbor in the case of resampling method = "near", or overshoots in the case of "cubic" methods.).This tool is intended to preserve full range of values when this is deemed more important.
#' @export
#'
#' @examples data("denali")
#' denali <- denali
#' denali <- toraster(denali)
#' minmax(denali)
#' #reproject to a coarser resolution
#' denali0 <- reproject(denali, lat = 63, lon = -151, rs = 1000,  h=100000, w=200000)
#' #note values less extreme than original
#' minmax(denali0)
#' denali1 <- RestoreMaxMin(denali, denali0)
#' #note max/min values are restored
#' minmax(denali1)
#'
#' #supplemental point elevations
#' data("peaks")
#' peaks <- peaks
#' #reproject to a finer resolution
#' denali0 <- reproject(denali, lat = 63, lon = -151, rs = 100,  h=100000, w=200000)
#' peaks <- st_transform(peaks, crs(denali0))
#' plot(denali0); plot(st_geometry(peaks), add=TRUE)
#' #Note that highest peak is known to be higher than raster indicates even if preserving original raster values when projected at higher resolution.
#' peaks[peaks$name %in% 'Denali',]$summit
#' minmax(denali0)
#' denali1 <- RestoreMaxMin(x=denali, y=denali0, s=peaks, e='summit')
#' #New high point is equal to supplemental point data.
#' minmax(denali1)
RestoreMaxMin <- function(x,y,s=NA,e=NA){#experimental version that incorporates point data and maintains neighborhood relative to original resolution
  hfactor1 <- ifelse(terra::linearUnits(x) == 0, 111111.1, terra::linearUnits(x))
  hfactor2 <- ifelse(terra::linearUnits(y) == 0, 111111.1, terra::linearUnits(y))
  rs <- pmax(res(x)[1]*hfactor1, res(y)[1]*hfactor2)
  #dumb raster
  y.rast <- rast(xmin=ext(y)[1]-rs, xmax=ext(y)[2]+rs,
                 ymin=ext(y)[3]-rs, ymax=ext(y)[4]+rs, crs=crs(y), res=rs)
  xmax <- project(x, y.rast, method='max')
  xmin <- project(x, y.rast, method='min')
  if(!is.na(e)){#incorporate point dataset if available
    pts <- sf::st_transform(s, crs=crs(y))
    pts <- rasterize(pts, y.rast, field = e)
    xmax <- max(xmax, pts, na.rm = TRUE)
  }

  xmax <- focalmax(xmax, r=rs*2, p='medium')
  xmin <- focalmin(xmin, r=rs*2, p='medium')

  #extremes to allow in target raster based on original raster
  supermax <- focalmax(xmax, r=rs*3, p='medium')
  supermin <- focalmin(xmin, r=rs*3, p='medium')
  #bring layers to target projection
  supermax <- project(supermax, y, method='bilinear')
  supermin <- project(supermin, y, method='bilinear')
  xmax <- project(xmax, y, method='bilinear')
  xmin <- project(xmin, y, method='bilinear')


  #wider neighborhood of target raster to avoid edge artifacts when using neighborhood of original raster
  ymean <- focalmed(y, r=rs*3, p='medium')
  ymax <- focalmax(y, r=rs*3, p='medium')
  ymin <- focalmin(y, r=rs*3, p='medium')
  #weave in maxmin of target raster.
  xmax <- max(xmax, ymax, na.rm = TRUE)
  xmin <- min(xmin, ymin, na.rm = TRUE)
  #push back any overshoots from target raster
  xmax <- min(xmax, supermax, na.rm = TRUE)
  xmin <- max(xmin, supermin, na.rm = TRUE)

  z <- ifel(ymax-ymean == 0 | ymean-ymin == 0, y,
            ifel(y >= ymean,
                 (y - ymean)/(ymax - ymean)*(xmax - ymean)+ymean,
                 (y - ymean)/(ymean - ymin)*(ymean - xmin)+ymean))

  return(z)
}








#' Create hillshade directly from DEM raster
#'
#' @param x Digital elevation model
#' @param angle The elevation angle(s) of the light source (sun), in degrees
#' @param direction The direction (azimuth) angle(s) of the light source (sun), in degrees
#'
#' @return Hillshade raster for use to visualize terrain texture.
#' This function combines terra::terrain() function to create slope and aspect rasters, then employs the terra::shade() function to generate hillshade in one step.
#' @export
#'
#' @examples data(denali)
#' denali <- denali
#' denali <- toraster(denali)
#' hsd <- hillshade(denali)
#' plot(hsd)
hillshade <- function(x, angle=45, direction=0){
  require(terra)
  hsd <-  shade(slope=terrain(x, 'slope', unit="radians"), aspect= terrain(x, 'aspect', unit="radians"), angle=angle, direction = direction)
  return(hsd)
}


#' From raster to list object
#'
#' @param x Terra raster object.
#'
#' @return List object suitable for saving as package data.
#' @export
#'
#' @examples
fromraster <- function(x){
  require(terra)
  y <- list(data = as.matrix(x,wide=TRUE),
            crs = crs(x),
            ext = t(as.matrix(ext(x))),
            names = names(x))
  return(y)}

#' To raster from list object
#'
#' @param x List object imported from package data.
#'
#' @return Terra raster object.
#' @export
#'
#' @examples data(denali)
#' denali <- denali
#' denali <- toraster(denali)
#' plot(denali)
toraster <- function(x){
  require(terra)
  y <- rast(x$data, crs = x$crs, ext = ext(x$ext))
  names(y) <- x$names
  return(y)}

makeWedge <- function(x, r, asp){
  u <- linearUnits(x)
  u <- ifelse(u == 0, 10000000/90,u)
  rs <- res(x)[1]*u
  size <- floor(r/rs)*2+1
  center <- floor(size/2)+1
  m <- matrix(1, nc=size, nr=size)
  y = 2*(center-row(m))/size
  x = 2*(col(m)-center)/size
  d <- ((y-0)^2+(x-0)^2)^0.5
  a <- ifelse(x < 0, 360-(90-asin(y/d)*360/2/pi),(90-asin(y/d)*360/2/pi))
  s <- 1-d

  asp <- ifelse(asp > 180, asp - 360, asp)
  a <- ifelse(a > 180, a - 360, a)
  a2 <- ifelse(abs(asp - a) > 180, 360-abs(asp - a), abs(asp - a))
  s1 <- s*(45-a2)/45
  s1 <- ifelse(a2 > 45 | s < 0,0, s1)
  s1 <- ifelse(row(m) == center & col(m) == center,1,s1)
  s1 <- s1^0.5
  return(s1)
}

#' Directional Maximum
#'
#' @param x raster
#' @param r radius (meters)
#' @param asp aspect or direction of the wind.
#' @param p precision, from low to exact, with higher levels of precision requiring more processing time
#'
#' @returns Raster with directional maximum with effects tapering to zero at radius, spread across 45 degrees. Subtracting elevation from directional maximum results in the equivalent to a rain shadow in the direction opposite of the specified aspect or wind direction.
#' @export
#'
#' @examples
shadowMax <- function(x, r, asp=270, p=c('low', 'medium', 'high','exact')){
  require(terra)
  #establish aggregating factor when radius is too large
  #p is for precision options
  u <- linearUnits(x)
  u <- ifelse(u == 0, 10000000/90,u)
  rs <- res(x)[1]*u
  p=p[1]
  if(p == 'low'){
    fc = floor(r/rs/9+1)
  }else if(p == 'medium'){
    fc = floor(r/rs/21+1)
  }else if(p == 'high'){
    fc = floor(r/rs/81+1)
  }else{
    fc = floor(r/rs/243+1)}

  x1 = x
  if(fc > 1){
    x1 <- aggregate(x, fact = fc, fun = 'max',  na.rm=TRUE)
  }
  #create a focal weights matrix of the appropriate size for given radius
  #fm <- focalMat(x1, d=r, type = 'circle')
  fm <- makeWedge(x1, r=r, asp=asp)

  #reduce mat size if raster too small
  matsize = nrow(fm)
  centermat <- floor(nrow(fm)/2+1)
  goodrows <- intersect((1:matsize),(centermat-floor(nrow(x1)-1)):(centermat+floor(nrow(x1)-1)))
  goodcols <- intersect((1:matsize),(centermat-floor(ncol(x1)-1)):(centermat+floor(ncol(x1)-1)))
  fm <- fm[goodrows,goodcols, drop = FALSE]
  if(nrow(fm) <= 1 & ncol(fm) <= 1){
    x1.max <- x1
  }else{
    x1.max <- focal(x1, fm, fun='max', na.rm=T)}
  #restore resolution in result
  if(fc > 1){
    x1.max <- project(x1.max, x)
  }
  return(x1.max)}



#' Rotate XY Coordinates
#'
#' @param x vector of x coordinates
#' @param y vector of y coordinates
#' @param a angle in degrees
#' @param cx optional center of rotation x coordinate (default is center of point cloud)
#' @param cy optional center of rotation y coordinate (default is center of point cloud)
#'
#' @returns data frame of rotated xy coordinates
#' @export
#'
#' @examples df <- data.frame(
#' @examples x=runif(10,0,10),
#' @examples y=rnorm(10,5,5))
#' @examples df2 <-  rotatexy(df$x,df$y, a=2)
#' @examples plot(df$y ~ df$x)
#' @examples points(df2$y ~ df2$x, col='red')
rotatexy <- function(x, y, a, cx = NA, cy = NA){
  df <- data.frame(x=x,y=y)

  if(is.na(cx) | is.na(cy)){
    cx <- mean(df$x)
    cy <- mean(df$y)}

  df$y0 <- df$y-cy
  df$x0 <- df$x-cx
  df$h <- ((df$x0)^2+(df$y0)^2)^0.5
  df$a0 <- ifelse(df$h==0,0,acos(df$y0/df$h))
  a1 <- a/360*2*pi
  df$a0 <- ifelse(df$x0 >= 0,df$a0,-1*df$a0)
  xr = ifelse(df$h==0,0,df$h*sin(df$a0+a1))+cx
  yr = ifelse(df$h==0,0,df$h*cos(df$a0+a1))+cy

  rdf <- data.frame(x=xr,y=yr)
  return(rdf)
}


#' Make XY Raster
#'
#' @param x raster to extract xy coordinates
#' @param rotations Specify number of rotations of xy coordinates (default zero)
#'
#' @returns Multi channel raster with xy coordinates. (rotations add alternative angles a random forest covariates). Rasters named "lat" for y coordinate, and "lon"  for x coordinate.
#' @export
#'
#' @examples x <- terra::rast(matrix(1:25, nrow=5, ncol=5))
#' @examples xyrast <- makexyrast(x)
makexyrast <- function(x, rotations=0){
  require(terra)
  angles <- 90/(rotations+1)
  df <- terra::as.data.frame(x, xy=TRUE)
  lat <- terra::rast(df[,c("x","y","y")],type="xyz", crs=crs(x)); names(lat)='lat'
  lon <- terra::rast(df[,c("x","y","x")],type="xyz", crs=crs(x)); names(lon)='lon'
  xyrast <- c(lat,lon)
  if(rotations > 0){
    for(i in 1:rotations){
      df0 <- climatools::rotatexy(df$x,df$y, a=angles*i)
      lat0 <- terra::rast(cbind(df[,c("x","y")], df0$y),type="xyz", crs=crs(x))
      lon0 <- terra::rast(cbind(df[,c("x","y")], df0$x),type="xyz", crs=crs(x))
      namelat <- paste0("lat",i)
      namelon <- paste0("lon",i)
      names(lat0) <- namelat
      names(lon0) <- namelon
      assign(namelat,lat0)
      assign(namelon,lon0)
    }
    xyrast <- terra::rast(mget(c("lat","lon",paste0("lat",1:rotations), paste0("lon",1:rotations))))
  }

  return(xyrast)
}







#' Re-fit climate later to DEM
#'
#' This function can be used to omit problematic portions of a climatic layer and refit to a digital elevation model. The results can be used to patch problems in larger climatic layers. Omitting elevation layer bases refit exclusively on trends in xy coordinates.
#'
#' @param x climatic layer to re-fit.
#' @param elev optional digital elevation model of similar resolution as climatic layer.
#' @param cropto spatial extent of cropped area (xmin, xmax, ymin, ymax)
#' @param cropfrom spatial extent of problematic area to be removed from model (xmin, xmax, ymin, ymax)
#' @param rotations Optional number of rotation of XY coordinates used for smoother random forest model of linear model residuals.
#' @param sampdens sample density for extracting points to train models.
#' @param altlayer optional alternative covariate layer(s) (can be multiple rasters combined, e.g. water body influence layer, aspects)
#'
#' @returns Re-fitted climatic model.
#' @export
#'
#' @examples
refit <- function(x, elev=NULL, cropto=NULL, cropfrom=NULL, rotations=0, sampdens = 1500, altlayer = NULL){
  require(terra)

  t0 <- x[[1]]
  #If cropping extent provided
  if(!is.null(cropto)){
    t0 <- crop(x[[1]], cropto);
  }

  #determine units of layer and rescale focal analyses so that they are proportional to resolution
  u <- terra::linearUnits(t0)
  u <- ifelse(u == 0, 10000000/90, u)
  rs <- (terra::res(t0)*u)[1]

  #Standardize name for extracting to a data frame to be used by formula
  # names(t0) <- "t0"

  #dummy values for having no data for these layers
  erel <- NULL
  wt1 <- NULL; wt2 <- NULL; wt3 <- NULL

  #Create alternative rotated XY coordinates to make a smoother random forest model.
  xy0 <- makexyrast(t0,rotations)

  #Omit problematic data to patch with model using less problematic data from adjacent area.
  if(!is.null(cropfrom)){
    t00 <- crop(x, cropfrom)
    t00 <- ifel(t00>0,NA,NA)
    t0 <- terra::merge(t00,t0, na.rm=FALSE)
  }

  #If elevation data is used, create auxiliary layer for relative elevation which captures local inversions.
  if(!is.null(elev)){
    # names(elev) <- 'elev'
    #get elevation data to match extent and resolution of input data
    e0 <- project(elev, t0) |> crop(ext(t0))
    # names(e0) <- 'elev'
    #generate relative elevation model
    emd <- focalmed(e0, rs*10)
    erel <- e0-emd
    erel <- ifel(erel > 0,erel,0)^0.5 - ifel(-erel > 0,-erel,0)^0.5; names(erel)<-'erel'
    rss <- c(t0, e0, erel, xy0)
    df0 <- terra::spatSample(c(t0, e0, erel, xy0), size=sampdens, xy=FALSE, values=TRUE)
    wt1 <- 1; wt2 <- 1
  }else{
    #Use only xy data instead
    rss <- c(t0, xy0)
    df0 <- terra::spatSample(c(t0, xy0), size=sampdens, xy=FALSE, values=TRUE)
  }
  if(!is.null(altlayer)){
    altlayer <- project(altlayer, t0) |> crop(ext(t0))
    # names(altlayer) <- 'altlayer'
    wt3 <- (1:length(names(altlayer)))*0+1
    rss <- c(rss, altlayer)
    df0 <- terra::spatSample(c(t0, e0, erel, xy0, altlayer), size=sampdens, xy=FALSE, values=TRUE)
  }

  #Standardize name for extracting to a data frame to be used by formula
  # names(t0) <- "t0"

  #Define formulas for general linear and random forest models.
  depvar <- names(t0)
  covars1 <- c(names(xy0)[1:2],names(elev),names(erel),names(altlayer))
  covars2 <- c(names(xy0),names(elev),names(erel),names(altlayer))
  wts <- c(c(1:((rotations+1)*2))*0+1/(rotations+1), wt1,wt2,wt3)
  covars3 <- c(covars1,'resids')

  #Sample rasters to train models (omitting missing data)
  # df0 <- terra::spatSample(c(t0, e0, erel, xy0), size=sampdens, xy=TRUE, values=TRUE)
  df0 <- df0[!is.na(df0[,depvar]),]

  #data points converted to spatial features to test coverage
  #dfsp <- sf::st_as_sf(df0, coords = c(x='lon', y='lat'), crs=sf::st_crs(t1)); plot(vect(dfsp))

  f.glm <- stats::as.formula(paste(paste(depvar,paste(paste(covars1, collapse = " + ", sep = ""),""), sep = " ~ ")
  ))

# environment(f.glm) <- environment()

  f.rf <- stats::as.formula(paste(paste("resids",paste(paste(covars2, collapse = " + ", sep = ""),""), sep = " ~ ")
  ))

  f.glm2 <- stats::as.formula(paste(paste(depvar,paste(paste(covars3, collapse = " + ", sep = ""),""), sep = " ~ ")
  ))

  #Run initial model to get linear trends
  gm <- stats::glm(f.glm,
            family='gaussian',
            data=df0)

    df0$pred <- predict(gm, df0)
    df0$resids <- df0[,depvar]-df0$pred

  # #Create additional layer with residuals using random forest model.
  rf <- ranger::ranger(f.rf,
               split.select.weights=wts,
               #num.trees = 1500,
               data=df0)
  resids <- terra::predict(rss, rf)
  resids <- focalmed(resids, rs*3); names(resids) <- 'resids'
  rss1 <- c(rss, resids)

if(is.null(elev) & is.null(altlayer)){
  df1 <- terra::spatSample(c(t0, xy0,resids), sampdens)
}else if(is.null(altlayer)){
  df1 <- terra::spatSample(c(t0, e0, erel, xy0,resids), sampdens)
}else{
  df1 <- terra::spatSample(c(t0, e0, erel, xy0, altlayer,resids), sampdens)
}
  df1 <- df1[!is.na(df1[,depvar]),]

  gm2 <- glm(f.glm2,
             family='gaussian',
             data=df1)

  pred <- terra::predict(rss1, gm2) ; names(pred) <- 'pred'
  return(pred)

  #nameing within a raster enclosed within a package function doesn't always work for spatSample
  #spatSample embedded within a package function only works directly with raster objects and can only concatenate within the function.

  }













#' Make climate raster from point data
#'
#' @param pts data frame of 3 columns; first two: xy coordinates (longitude and latitude decimal degrees); third column is climate attribute like temperature.
#' @param altlayer terra raster layer (multiple layers).
#' @param cropto Optional crop extent c(xmin, xmax, ymin, ymax).
#' @param covrange Minimum range in first covariate before building submodel coefficients.
#' @param minrow Minimum number of data points before building submodel coefficients.
#' @param segx Number of segments in x axis for building submodel coefficients.
#' @param segy Number of segments in y axis for building submodel coefficients.
#' @param cropbuffer Buffer around the crop extent to incorporate points into model outside cropped area (units of raster projection).
#' @param randforest Use random forest model if true (requires ranger package) to generate coefficients and residuals layers. Alternatively, coefficients and residuals layers are interpolated solely using xy coordinate using gstats package.
#'
#' @returns Climate raster matching cropped extent.
#' @export
#'
#' @examples
toclimrast <- function(pts, altlayer, cropto=NULL, covrange=0, minrow=50, segx=5, segy=5,
                       cropbuffer=5, randforest = TRUE){
  #crop full extent if null
  if(is.null(cropto)){cropto <- terra::ext(altlayer)
  cropbuffer=0}
  #ensure that input has only 3 columns
  pts <- as.data.frame(pts)
  pts <- pts[,1:3]
  names(pts) <- c('x','y','z')

  #crop raster to new extent with buffer
  cropto0 <- cropto + c(-cropbuffer,cropbuffer,-cropbuffer,cropbuffer)
  #crop point data extent and convert to terra spatial vector.
  pts <- pts |> subset(x >= cropto0[1] & x <= cropto0[2] &
                         y >= cropto0[3] & y <= cropto0[4])
  vts <- vect(pts, geom=c("x", "y"),crs=crs('epsg:4326'))

  grd <- crop(altlayer, cropto0)
  #add raster with rotated xy coordinates
  xy0 <- climatools::makexyrast(grd[[1]],6)
  grdall <- c(grd, xy0)
  #create lower resolution raster to model covariates and residuals
  grdall.1 <- aggregate(grdall,fact=3,fun="mean")
  #get raster units to ensure that focal analyses neighborhoods are consistent
  u <- terra::linearUnits(grd)
  u <- ifelse(u == 0, 10000000/90, u)
  rs <- (terra::res(grd)*u)[1]

  #extract rasters to points
  vts <- project(vts,altlayer)
  vtsgrd <- terra::extract(grdall,vts)
  pts <- cbind(pts,vtsgrd)

  #build formulas
  depvar <- names(pts)[3]
  covars1 <- c(names(grd),names(xy0)[1:2])
  covars2 <- c(names(grd),names(xy0))
  f.glm <- stats::as.formula(paste(paste(depvar,paste(paste(covars1, collapse = " + ", sep = ""),""), sep = " ~ ")
  ))

  #prepare regression loops
  exfactors <- c(0,0.5,1,2,5,10)
  spanx <- (cropto[2]-cropto[1])/segx
  spany <- (cropto[4]-cropto[3])/segy
  pts$inner <- NA
  pts$outer <- NA
  pts$coeffs0 <- NA
  pts$erange <- NA
  #loop through x and y segments
  for(i.x in 1:segx){
    for(i.y in 1:segy){
      #set status for this segment until acceptable regression model
      success <- FALSE
      #gradually expand size of analysis area
      for(i.f in 1:length(exfactors)){
        #i.x=3;i.y=4;i.f=1
        if(!success){
          exfact <- exfactors[i.f]
          addtoy <- spany*exfact
          addtox <- spanx*exfact
          #establish inner and outer points; inner points of segment carry the values of covariates; outer values are involved in the building the models
          crop0 <- c(cropto[1]+(i.x-1)*spanx,cropto[1]+i.x*spanx,
                     cropto[3]+(i.y-1)*spany,cropto[3]+i.y*spany)
          pts <- pts |> mutate(inner = ifelse(x >= crop0[1] & x <= crop0[2] &
                                                y >= crop0[3] & y <= crop0[4], 1, 0),
                               outer = ifelse(x >= (crop0[1]-addtox) & x <= (crop0[2]+addtox) &
                                                y >= (crop0[3]-addtoy) & y <= (crop0[4]+addtoy), 1, 0))
          pts.i <- pts |> subset(outer ==1)
          #move on if segment is has too few rows to build model
          if(nrow(pts.i) > minrow){
            #move on (expand extent) if segment points do not have sufficient range in first variable (e.g. elevation) to build accurate model
            erange0 <- max(pts.i[,5])-min(pts.i[,5])
            if(erange0 >= covrange){
              #model segment of points and feed coefficients into points dataset
              gm <- stats::glm(f.glm,
                               family='gaussian',
                               data=pts.i)
              summary(gm)
              cofs <- list(gm$coefficients)
              pts <- pts |> mutate(coeffs0 = ifelse(inner %in% 1, cofs,coeffs0),
                                   erange = ifelse(inner %in% 1, erange0,erange))
              success <- TRUE}
          }}}
    }}

  #extract coefficients from points and convert to rasters using either randomforest model or interpolation
  cflist <-t(as.data.frame(pts$coeffs0))
  nc <- ncol(cflist)
  grdall.0 <- rast(resolution=res(grdall.1), crs=crs(grdall.1), extent=ext(grdall.1), nlyrs=nc)
  for(i in 1:nc){
    pts$coeffs <- cflist[,i]
    if(randforest){
      f.rf <- stats::as.formula(paste(paste("coeffs",paste(paste(covars2, collapse = " + ", sep = ""),""), sep = " ~ ")
      ))
      rf <- ranger::ranger(f.rf,
                           # split.select.weights=wts,
                           #num.trees = 1500,
                           data=pts[!is.na(pts$coeffs),])

      cofffs <- terra::predict(grdall.1, rf)
    }else{
      xyz <- pts[,c('x','y','coeffs')] |> subset(!is.na(coeffs))
      gs <- gstat::gstat(formula=coeffs~1, locations=~x+y, data=xyz, nmax=32, set=list(idp = 2))
      cofffs <- interpolate(grdall.1, gs, debug.level=0)[[1]]
    }
    cofffs <- focalmed(cofffs, segy*u/3);
    names(cofffs)  <- paste0("coef.",i)
    grdall.0[[i]] <- cofffs
  }
  grdall.0 <- project(grdall.0, grdall)

  #extract coefficients to points
  pts2 <- pts |> cbind(extract(grdall.0,vts))
  grdall2 <- c(grdall, grdall.0)
  #create formula with covariates and coefficients
  covarc1 <- names(grdall.0)[2:nc]
  intcp <- names(grdall.0)[1]
  f.glm2 <- stats::as.formula(paste(depvar,paste(intcp, paste(covars1,"*",covarc1, collapse = " + ", sep = ""), sep = " + "), sep = " ~ "))

  #linear model with new formula
  gm2 <- stats::glm(f.glm2,
                    family='gaussian',
                    data=pts2)
  summary(gm2)
  1-gm2$deviance/gm2$null.deviance
  #use model to generate prediction layer
  pred <- terra::predict(grdall2, gm2)
  #use model to generate residuals in points
  pts2$pred <- predict(gm2,pts2)
  pts2$resid <- pts2$z-pts2$pred
  #build formula to generate residual raster with either randomforest model or interolation
  if(randforest){
    f.rf2 <- stats::as.formula(paste(paste("resid",paste(paste(covars2, collapse = " + ", sep = ""),""), sep = " ~ ")
    ))
    rf2 <- ranger::ranger(f.rf2,
                          # split.select.weights=wts,
                          #num.trees = 1500,
                          data=pts2[!is.na(pts2$resid),])
    resid <- terra::predict(grdall.1, rf2)
  }else{
    xyz <- pts2[,c('x','y','resid')]
    gs <- gstat::gstat(formula=resid~1, locations=~x+y, data=xyz, nmax=32, set=list(idp = 2))
    resid <- interpolate(grdall.1, gs, debug.level=0)[[1]]
  }

  #add residual layer to linear model prediction layer
  resid <- resid |> climatools::focalmed(50000)  |> project(grdall)
  model <- terra::crop(resid+pred, cropto)
  return(model)
}













