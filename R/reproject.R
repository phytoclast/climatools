# library(climatools)
# library(terra)
# library(sf)
# dem <- rast('C:/a/geo/dem/SRTM250m/usa/w001001.adf')
# dem <- rast('C:/workspace2/lidar/output/mio1/canopy.tif')

#reproject raster to metric units centered to middle of raster ----

#' Set a custom projection.
#'
#' @param prj Projection method among 12 options (default is equal area azimuthal).
#' @param datum Set of three common datums (default is WGS84).
#' @param unit Set units of projection (default is meters)
#' @param lat Latitude of center of projection, first parallel of a conic projection, or first point of an oblique Mercator projection.
#' @param lon Longitude of center of projection, center meridian of a conic projection (for scaling), or first point of an oblique Mercator projection.
#' @param lat2 Latitude of second parallel of a conic projection (for scaling), or second point of an oblique Mercator projection (to establish scaling axis).
#' @param lon2 Longitude of second point of an oblique Mercator projection (to establish scaling axis).
#' @param feast False easting.
#' @param fnorth False northing.
#' @param orglat Latitude of coordinate origins.
#' @param scf Scaling factor of Transverse and Oblique Mercator projections.
#'
#' @return Well Known Text projection string to be used in raster and vector coordinate reference definitions.
#' @export
#'
#' @examples
setProjection <- function(prj = c('point.equalarea','point.equaldistant','conformal.point',
                                  'cylinder.equalarea','cylinder.equaldistant','conformal.cylindric',
                                  'cone.equalarea','cone.equaldistant','conformal.conic',
                                  'conformal.transverse','conformal.oblique','geographic'),
                          datum = c('WGS84','NAD83','NAD27'),unit=c('m','ft','usft'),
                          lat=NA, lon=NA, lat2=NA,lon2=NA,feast=0,fnorth=0, orglat=NA, scf=0.9996){
  prj <- prj[1];datum <- datum[1];unit <- unit[1]
  #alternative numeric parameter assigned to text
  if(is.numeric(prj)){
    preprj <- c('point.equalarea','point.equaldistant','conformal.point',
                'cylinder.equalarea','cylinder.equaldistant','conformal.cylindric',
                'cone.equalarea','cone.equaldistant','conformal.conic',
                'conformal.transverse','conformal.oblique','geographic')
    prj <- preprj[prj]
  }
  if(!is.numeric(datum)){
    predatum <- c('WGS84','NAD83','NAD27')
    datum <- which(predatum %in% datum)
  }
  if(!is.numeric(unit)){
    preunit <- c('m','ft','usft')
    unit <- which(preunit %in% unit)
  }

  da <- datum

  un <- unit


  begin <- 'PROJCS["Custom",'

  datums <- c(
    #WGS 84
    'GEOGCS["WGS 84",
        DATUM["WGS_1984",
            SPHEROID["WGS 84",6378137,298.257223563]],
        PRIMEM["Greenwich",0,
            ANGLEUNIT["degree",0.0174532925199433]],
        ID["EPSG",4326]],',
    #NAD83
    'BASEGEOGCRS["NAD83",
        DATUM["North American Datum 1983",
            ELLIPSOID["GRS 1980",6378137,298.257222101,
                LENGTHUNIT["metre",1]]],
        PRIMEM["Greenwich",0,
            ANGLEUNIT["degree",0.0174532925199433]],
        ID["EPSG",4269]],',
    #NAD27
    'BASEGEOGCRS["NAD27",
            DATUM["North American Datum 1927",
                  ELLIPSOID["Clarke 1866",6378206.4,294.978698213898,
                            LENGTHUNIT["metre",1]]],
            PRIMEM["Greenwich",0,
                   ANGLEUNIT["degree",0.0174532925199433]],
            ID["EPSG",4267]],'
  )
  datums2 <- c('EPSG:4326','EPSG:4269','EPSG:4267')


  unitss <- c('UNIT["metre",1,AUTHORITY["EPSG","9001"]]',
              'UNIT["foot",0.3048, AUTHORITY["EPSG","9002"]]',
              'UNIT["US survey foot",0.304800609601219,AUTHORITY["EPSG","9003"]]')



  axises <- ',AXIS["Easting",EAST], AXIS["Northing",NORTH]]]'


  projections <- c(
    #1.azimutal 1,2
    'PROJECTION["Lambert_Azimuthal_Equal_Area"],',
    'PROJECTION["Azimuthal_Equidistant"],',
    #2.stereographic 3
    'PROJECTION["Stereographic"],',
    #3.cylindric 4,5
    'PROJECTION["Cylindrical_Equal_Area"],',
    'PROJECTION["Equirectangular"],',
    #4.mercator 6
    'PROJECTION["Mercator_2SP"],',
    #5.conic 7,8
    'PROJECTION["Albers_Conic_Equal_Area"],',
    'PROJECTION["Equidistant_Conic"],',
    #6.conformal conic 9
    'PROJECTION["Lambert_Conformal_Conic_2SP"],',
    #7.transverse 10
    'PROJECTION["Transverse_Mercator"],',
    #8.oblique 11
    'PROJECTION["Hotine_Oblique_Mercator_Two_Point_Natural_Origin"],'
  )



  if(is.na(lat2)){lat1 <- lat - 8}else{lat1 <- lat}
  if(is.na(lon2)){lon1 <- lon - 3}else{lon1 <- lon}
  if(is.na(lat2)){lat2 <- lat + 8}
  if(is.na(lon2)){lon2 <- lon + 3}
  if(is.na(orglat)){orglat <- lat}
  lat <- mean(lat1,lat2);lon <- mean(lon1,lon2)

  azcentlat <- paste0('PARAMETER["latitude_of_center",',lat,'],')
  azcentlon <- paste0('PARAMETER["longitude_of_center",',lon,'],')
  parfeast <- paste0('PARAMETER["false_easting",',feast,'],')
  parfnorth <- paste0('PARAMETER["false_northing",',fnorth,'],')
  parscf <- paste0('PARAMETER["scale_factor",',scf,'],')


  coniclat1 <- paste0('PARAMETER["standard_parallel_1",',lat1,'],')
  coniclat2 <- paste0('PARAMETER["standard_parallel_2",',lat2,'],')
  coniclat0 <- paste0('PARAMETER["latitude_of_origin",',orglat,'],')
  coniclon <- paste0('PARAMETER["central_meridian",',lon,'],')
  coniclatcenter <- paste0('PARAMETER["latitude_of_center",',lat,'],')
  conicloncenter <- paste0('PARAMETER["longitude_of_center",',lon,'],')

  oblilat1 <- paste0('PARAMETER["latitude_of_point_1",',lat1,'],')
  oblilon1 <- paste0('PARAMETER["longitude_of_point_1",',lon1,'],')
  oblilat2 <- paste0('PARAMETER["latitude_of_point_2",',lat2,'],')
  oblilon2 <- paste0('PARAMETER["longitude_of_point_2",',lon2,'],')


  ifaz <- paste0(azcentlat, azcentlon)
  ifcyl <- paste0(coniclat1, coniclon)
  ifcon <- paste0(coniclatcenter, conicloncenter, coniclat1, coniclat2)
  ifconf <- paste0(coniclat0, coniclon, coniclat1, coniclat2)
  ifutm <- paste0(coniclat0, coniclon, parscf)
  ifobl <- paste0(azcentlat, oblilat1, oblilon1, oblilat2, oblilon2, parscf)
  ifster <- paste0(coniclat0, coniclon)
  ifmerc <- paste0(coniclat1,  coniclon)

  pickparameters <- c(ifaz, ifcyl, ifcon, ifconf, ifutm, ifobl, ifster, ifmerc)



  if(prj %in% 'point.equalarea'){
    prochoice <- 1
    pramchoice <- 1
  }
  if(prj %in% 'point.equaldistant'){
    prochoice <- 2
    pramchoice <- 1
  }
  if(prj %in% 'conformal.point'){
    prochoice <- 3
    pramchoice <- 7
  }

  if(prj %in% 'cylinder.equalarea'){
    prochoice <- 4
    pramchoice <- 2
  }
  if(prj %in% 'cylinder.equaldistant'){
    prochoice <- 5
    pramchoice <- 2
  }
  if(prj %in% 'conformal.cylindric'){
    prochoice <- 6
    pramchoice <- 2
  }
  if(prj %in% 'cone.equalarea'){
    prochoice <- 7
    pramchoice <- 3
  }
  if(prj %in% 'cone.equaldistant'){
    prochoice <- 8
    pramchoice <- 3
  }
  if(prj %in% 'conformal.conic'){
    prochoice <- 9
    pramchoice <- 4
  }

  if(prj %in% 'conformal.transverse'){
    prochoice <- 10
    pramchoice <- 5
  }
  if(prj %in% 'conformal.oblique'){
    prochoice <- 11
    pramchoice <- 6
  }
  if(prj %in% 'geographic'){
    assembly <- crs(datums2[da])
  }else{
    assembly <- paste0(begin,
                       datums[da],
                       projections[prochoice],
                       pickparameters[pramchoice],
                       parfeast, parfnorth, unitss[un],axises)
  }

  return(assembly)

}


#' Reproject any raster to Lambert_Azimuthal_Equal_Area with specified center and extent.
#'
#' @param dem Digital elevation model or other type of raster
#' @param lat center latitude (default is existing raster)
#' @param lon center longitude (default is existing raster)
#' @param rs resolution of raster in meters (default is converted from existing raster)
#' @param h extent north to south in meters (default is converted from existing raster)
#' @param w extent east to west in meters (default is converted from existing raster)
#' @param prj Coordinate reference system (crs); Defaults to equal area azimuth.
#'
#' @return new raster with Lambert_Azimuthal_Equal_Area projection
#' @export
#'
#' @examples
reproject <- function(dem, lat=NA, lon=NA, rs = NA,  h = NA, w = NA, prj=NA){
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
if(is.na(prj)){
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
    UNIT["metre",1]]')}else{wkt.new <- prj}

  #find units and resolution of current projection
  ishfeet <- (grepl('unit.*"foot".*vertcrs',tolower(st_crs(dem))) | !grepl('vertcrs',tolower(st_crs(dem))) & grepl('unit.*"foot"',tolower(st_crs(dem))))
  ishusfeet <- (grepl('unit.*"us survey foot".*vertcrs',tolower(st_crs(dem)))| !grepl('vertcrs',tolower(st_crs(dem))) & grepl('unit.*"us survey foot"',tolower(dem)))
  isprojected <- grepl('conversion',tolower(st_crs(dem)))

  hfactor <- ifelse(!isprojected[length(st_crs(dem))], 111111.1,
                    ifelse(ishfeet[length(st_crs(dem))], 0.3048,
                           ifelse(ishusfeet[length(st_crs(dem))], 0.304800609601219, 1)))
  r = mean(res(dem))*hfactor


  #find extent after projection
  ex <- data.frame(rname=c('center','ul','ll','ur','lr'),
                   xcoord=c(lon, ext(dem)[1],ext(dem)[1],ext(dem)[2],ext(dem)[2]),
                   ycoord=c(lat, ext(dem)[3],ext(dem)[4],ext(dem)[3],ext(dem)[4]))

  ex <- sf::st_as_sf(as.data.frame(ex), coords = c("xcoord","ycoord"), crs=st_crs(dem))
  ex.trans <- (st_transform(ex,crs=wkt.new))
  ex.trans <- as.data.frame(st_coordinates(ex.trans))

  #set final extent
  if(is.na(h)){
    ymn= min(ex.trans$Y)
    ymx= max(ex.trans$Y)
  }else{
    ymn= ex.trans$Y[1]-h/2
    ymx= ex.trans$Y[1]+h/2
  }
  if(is.na(w)){
    xmn= min(ex.trans$X)
    xmx= max(ex.trans$X)
  }else{
    xmn= ex.trans$X[1]-w/2
    xmx= ex.trans$X[1]+w/2
  }

  #set final resolution
  r <- ifelse(is.na(rs),r,rs)

  #dumb raster
  y.rast <- rast(xmin=xmn, xmax=xmx,
                 ymin=ymn, ymax=ymx, crs=wkt.new, res=r)
  #project raster
  dem2 <- project(dem, y.rast, method = 'bilinear')

  return(dem2)}
