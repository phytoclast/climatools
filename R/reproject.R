# library(climatools)
# library(terra)
# library(sf)
# dem <- rast('C:/a/geo/dem/SRTM250m/usa/w001001.adf')
# dem <- rast('C:/workspace2/lidar/output/mio1/canopy.tif')

#reproject raster to metric units centered to middle of raster ----

#' Reproject any raster to Lambert_Azimuthal_Equal_Area with specified center and extent.
#'
#' @param dem Digital elevation model or other type of raster
#' @param lat center latitude (default is existing raster)
#' @param lon center longitude (default is existing raster)
#' @param rs resolution of raster in meters (default is converted from existing raster)
#' @param h extent north to south in meters (default is converted from existing raster)
#' @param w extent east to west in meters (default is converted from existing raster)
#'
#' @return new raster with Lambert_Azimuthal_Equal_Area projection
#' @export
#'
#' @examples
reproject <- function(dem=dem, lat=NA, lon=NA, rs = NA,  h = NA, w = NA){
  require(terra)
  require(sf)

  #get lat/lon for center of dem by projecting to 'EPSG:4326'
  xcoord =(ext(dem)[1] + ext(dem)[2])/2
  ycoord =(ext(dem)[3] + ext(dem)[4])/2

  pt = data.frame(
    xcoord,ycoord
  )
  pt <- sf::st_as_sf(as.data.frame(pt), coords = c("xcoord","ycoord"), crs=st_crs(dem))
  pt.trans <- (st_transform(pt,crs=st_crs('EPSG:4326')))
  pt.trans <- st_coordinates(pt.trans)


  #set final lat/lon of projection
  if(is.na(lat)){
    lat <- pt.trans[,2]
  }
  if(is.na(lon)){
    lon <- pt.trans[,1]
  }

  wkt.new <- paste0('PROJCS["Centered Equal Area",
    GEOGCS["WGS 84",
        DATUM["WGS_1984",
            SPHEROID["WGS 84",6378137,298.257223563]],
        PRIMEM["Greenwich",0],
        UNIT["Degree",0.0174532925199433]],
    PROJECTION["Lambert_Azimuthal_Equal_Area"],
    PARAMETER["latitude_of_center",',lat,'],
    PARAMETER["longitude_of_center",',lon,'],
    PARAMETER["false_easting",0],
    PARAMETER["false_northing",0],
    UNIT["metre",1]]')

  #find units and resolution of current projection
  ishfeet <- (grepl('unit.*"foot".*vertcrs',tolower(st_crs(dem))) | !grepl('vertcrs',tolower(st_crs(dem))) & grepl('unit.*"foot"',tolower(st_crs(dem))))
  ishusfeet <- (grepl('unit.*"us survey foot".*vertcrs',tolower(st_crs(dem)))| !grepl('vertcrs',tolower(st_crs(dem))) & grepl('unit.*"us survey foot"',tolower(dem)))
  isprojected <- grepl('conversion',tolower(st_crs(dem)))

  hfactor <- ifelse(!isprojected[length(st_crs(dem))], 111111.1,
                    ifelse(ishfeet[length(st_crs(dem))], 0.3048,
                           ifelse(ishusfeet[length(st_crs(dem))], 0.304800609601219, 1)))
  r = mean(res(dem))*hfactor


  #find extent after projection
  ex <- data.frame(rname=c('ul','ll','ur','lr'),
                   xcoord=c(ext(dem)[1],ext(dem)[1],ext(dem)[2],ext(dem)[2]),
                   ycoord=c(ext(dem)[3],ext(dem)[4],ext(dem)[3],ext(dem)[4])
  )
  ex <- sf::st_as_sf(as.data.frame(ex), coords = c("xcoord","ycoord"), crs=st_crs(dem))
  ex.trans <- (st_transform(ex,crs=wkt.new))
  ex.trans <- as.data.frame(st_coordinates(ex.trans))

  #set final extent
  if(is.na(h)){
    ymn= min(ex.trans$Y)
    ymx= max(ex.trans$Y)
  }else{
    ymn= -h/2
    ymx= h/2
  }
  if(is.na(w)){
    xmn= min(ex.trans$X)
    xmx= max(ex.trans$X)
  }else{
    xmn= -w/2
    xmx= w/2
  }

  #set final resolution
  r <- ifelse(is.na(rs),r,rs)

  #dumb raster
  y.rast <- rast(xmin=xmn, xmax=xmx,
                 ymin=ymn, ymax=ymx, crs=wkt.new, res=r)
  #project raster
  dem2 <- project(dem, y.rast, method = 'bilinear')

  return(dem2)}
