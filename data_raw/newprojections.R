begin <- 'PROJCS["Custom",'

datums <- c(
  #WGS 84
  'GEOGCS["WGS 84",
        DATUM["WGS_1984",
            SPHEROID["WGS 84",6378137,298.257223563]],
        PRIMEM["Greenwich",0],
        UNIT["Degree",0.0174532925199433]],',
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


units <- c('UNIT["metre",1,AUTHORITY["EPSG","9001"]]',
           'UNIT["foot",0.3048, AUTHORITY["EPSG","9002"]]',
           'UNIT["US survey foot",0.304800609601219,AUTHORITY["EPSG","9003"]]')



axises <- ',AXIS["Easting",EAST], AXIS["Northing",NORTH]]]'

#survey feet
'UNIT["degree",0.0174532925199433,
     AUTHORITY["EPSG","9122"]]'


projections <- c(
  #azimutal 1,2
  'PROJECTION["Lambert_Azimuthal_Equal_Area"],',
  'PROJECTION["Azimuthal_Equidistant"],',
  #stereographic 3
  'PROJECTION["Stereographic"],',
  #cylindric 4,5
  'PROJECTION["Cylindrical_Equal_Area"],',
  'PROJECTION["Equirectangular"],',
  #mercator 6
  'PROJECTION["Mercator_2SP"],',
  #conic 7,8
  'PROJECTION["Albers_Conic_Equal_Area"],',
  'PROJECTION["Equidistant_Conic"],',
  #conformal conic 9
  'PROJECTION["Lambert_Conformal_Conic_2SP"],',
  #transverse 10
  'PROJECTION["Transverse_Mercator"],',
  #oblique 11
  'PROJECTION["Hotine_Oblique_Mercator_Two_Point_Natural_Origin"],'
)
lat <- 40
lon <- -85
lat1 <- NA
lon1 <- NA
lat2 <- lat + 8
lon2 <- lon + 3

shape <- c('point', 'cylinder', 'cone')
equal <- c('area', 'distance', 'angle')
angle <- c('standard', 'transverse', 'oblique')

scf <- 0.9996
feast <- 0
fnorth <- 0
orglat <- NA



if(is.na(lat1)){lat1 <- lat - 8}
if(is.na(lon1)){lon1 <- lon - 3}
if(is.na(lat2)){lat2 <- lat + 8}
if(is.na(lon2)){lon2 <- lon + 3}
if(is.na(orglat)){orglat <- lat}

  
azcentlat <- paste0('PARAMETER["latitude_of_center",',lat,'],')
azcentlon <- paste0('PARAMETER["longitude_of_center",',lon,'],')
parfeast <- paste0('PARAMETER["false_easting",',feast,'],')
parfnorth <- paste0('PARAMETER["false_northing",',fnorth,'],')
parscf <- paste0('PARAMETER["scale_factor",',scf,'],')


coniclat1 <- paste0('PARAMETER["standard_parallel_1",',lat1,'],')
coniclat2 <- paste0('PARAMETER["standard_parallel_2",',lat2,'],')
coniclat0 <- paste0('PARAMETER["latitude_of_origin",',orglat,'],')
coniclon <- paste0('PARAMETER["central_meridian",',lon,'],')

oblilat1 <- paste0('PARAMETER["latitude_of_point_1",',lat1,'],')
oblilon1 <- paste0('PARAMETER["longitude_of_point_1",',lon1,'],')
oblilat2 <- paste0('PARAMETER["latitude_of_point_2",',lat2,'],')
oblilon2 <- paste0('PARAMETER["longitude_of_point_2",',lon2,'],')


ifaz <- paste0(azcentlat, azcentlon)
ifcyl <- paste0(coniclat1, coniclon)
ifcon <- paste0(azcentlat, azcentlon, coniclat1, coniclat2)
ifconf <- paste0(coniclat0, coniclon, coniclat1, coniclat2)
ifutm <- paste0(coniclat0, coniclon, parscf)
ifobl <- paste0(azcentlat, oblilat1, oblilon1, oblilat2, oblilon2, parscf)
ifster <- paste0(coniclat0, coniclon)
pickparameters <- c(ifaz, ifcyl, ifcon, ifconf, ifutm, ifobl, ifster)

assembly <- paste0(begin,
               datums[1],
               projections[1],
               pickparameters[1],
               parfeast, parfnorth, units[1],axises)





#WGS 84
st_crs('EPSG:4326')
'GEOGCS["WGS 84",
        DATUM["WGS_1984",
            SPHEROID["WGS 84",6378137,298.257223563]],
        PRIMEM["Greenwich",0],
        UNIT["Degree",0.0174532925199433]],'
#NAD83
st_crs('EPSG:4269')
'BASEGEOGCRS["NAD83",
        DATUM["North American Datum 1983",
            ELLIPSOID["GRS 1980",6378137,298.257222101,
                LENGTHUNIT["metre",1]]],
        PRIMEM["Greenwich",0,
            ANGLEUNIT["degree",0.0174532925199433]],
        ID["EPSG",4269]],'
#NAD27
st_crs('EPSG:4267')
'BASEGEOGCRS["NAD27",
            DATUM["North American Datum 1927",
                  ELLIPSOID["Clarke 1866",6378206.4,294.978698213898,
                            LENGTHUNIT["metre",1]]],
            PRIMEM["Greenwich",0,
                   ANGLEUNIT["degree",0.0174532925199433]],
            ID["EPSG",4267]],'
st_crs('EPSG:3857')
st_crs('EPSG:3079')
st_crs('EPSG:6202')
st_crs('EPSG:2251')
st_crs('ESRI:102008')
st_crs('ESRI:102010')
st_crs('EPSG:10598')

#equalarea
#equidistant
#conformal (equal angle)

#azimuthal equalarea
'PROJECTION["Lambert_Azimuthal_Equal_Area"],
        PARAMETER["latitude_of_center",52],
        PARAMETER["longitude_of_center",10],
        PARAMETER["false_easting",4321000],
        PARAMETER["false_northing",3210000],
        UNIT["metre",1,
            AUTHORITY["EPSG","9001"]],
        AXIS["Easting",EAST],
        AXIS["Northing",NORTH]]'

#cylindric equalarea
'PROJECTION["Cylindrical_Equal_Area"],
PARAMETER["standard_parallel_1",0],
PARAMETER["central_meridian",0],
PARAMETER["false_easting",0],
PARAMETER["false_northing",0],
UNIT["metre",1,
     AUTHORITY["EPSG","9001"]],
AXIS["Easting",EAST],
AXIS["Northing",NORTH]]'

#conic equalarea
'PROJECTION["Albers_Conic_Equal_Area"],
    PARAMETER["latitude_of_center",23],
    PARAMETER["longitude_of_center",-96],
    PARAMETER["standard_parallel_1",29.5],
    PARAMETER["standard_parallel_2",45.5],
    PARAMETER["false_easting",0],
    PARAMETER["false_northing",0],
    UNIT["metre",1,
        AUTHORITY["EPSG","9001"]],
    AXIS["Easting",EAST],
    AXIS["Northing",NORTH]]'

#oblique equalarea


#azimuthal equidistant
'PROJECTION["Azimuthal_Equidistant"],
    PARAMETER["latitude_of_center",90],
    PARAMETER["longitude_of_center",0],
    PARAMETER["false_easting",0],
    PARAMETER["false_northing",0],
    UNIT["metre",1,
        AUTHORITY["EPSG","9001"]],
    AXIS["Easting",EAST],
    AXIS["Northing",NORTH]'

#cylindric equidistant
'PROJECTION["Equirectangular"],
PARAMETER["standard_parallel_1",0],
PARAMETER["central_meridian",0],
PARAMETER["false_easting",0],
PARAMETER["false_northing",0],
UNIT["metre",1,
     AUTHORITY["EPSG","9001"]],
AXIS["Easting",EAST],
AXIS["Northing",NORTH]]'

#conic equidistant
'PROJECTION["Equidistant_Conic"],
PARAMETER["latitude_of_center",40],
PARAMETER["longitude_of_center",-96],
PARAMETER["standard_parallel_1",20],
PARAMETER["standard_parallel_2",60],
PARAMETER["false_easting",0],
PARAMETER["false_northing",0],
UNIT["metre",1,
     AUTHORITY["EPSG","9001"]],
AXIS["Easting",EAST],
AXIS["Northing",NORTH]'

#oblique equidistant

#azimuthal conformal
'PROJECTION["Polar_Stereographic"],
PARAMETER["latitude_of_origin",90],
PARAMETER["central_meridian",0],
PARAMETER["scale_factor",1],
PARAMETER["false_easting",0],
PARAMETER["false_northing",0]'


#cylindric conformal
'PROJECTION["Mercator_2SP"],
PARAMETER["standard_parallel_1",0],
PARAMETER["central_meridian",0],
PARAMETER["false_easting",0],
PARAMETER["false_northing",0],
UNIT["metre",1,
     AUTHORITY["EPSG","9001"]],
AXIS["Easting",EAST],
AXIS["Northing",NORTH]]'

#conic conformal
'PROJECTION["Lambert_Conformal_Conic_2SP"],
PARAMETER["latitude_of_origin",40],
PARAMETER["central_meridian",-96],
PARAMETER["standard_parallel_1",20],
PARAMETER["standard_parallel_2",60],
PARAMETER["false_easting",0],
PARAMETER["false_northing",0],
UNIT["metre",1,
     AUTHORITY["EPSG","9001"]],
AXIS["Easting",EAST],
AXIS["Northing",NORTH]]'

#transverse conformal
'    PROJECTION["Transverse_Mercator"],
    PARAMETER["latitude_of_origin",0],
    PARAMETER["central_meridian",15],
    PARAMETER["scale_factor",0.9996],
    PARAMETER["false_easting",500000],
    PARAMETER["false_northing",0],
    UNIT["metre",1,
        AUTHORITY["EPSG","9001"]],
    AXIS["Easting",EAST],
    AXIS["Northing",NORTH]]'

#oblique conformal
'PROJECTION["Hotine_Oblique_Mercator_Two_Point_Natural_Origin"],
    PARAMETER["latitude_of_center",40],
    PARAMETER["latitude_of_point_1",0],
    PARAMETER["longitude_of_point_1",0],
    PARAMETER["latitude_of_point_2",60],
    PARAMETER["longitude_of_point_2",60],
    PARAMETER["scale_factor",1],
    PARAMETER["false_easting",0],
    PARAMETER["false_northing",0],
'





reproject <- function(dem, lat=NA, lon=NA, rs = NA,  h = NA, w = NA){
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