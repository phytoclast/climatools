begining <- 'PROJCS["Custom",'

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
  #azimutal
  'PROJECTION["Lambert_Azimuthal_Equal_Area"],',
  'PROJECTION["Azimuthal_Equidistant"],',
  'PROJECTION["Stereographic"],',
  #cylindric
  'PROJECTION["Cylindrical_Equal_Area"],',
  'PROJECTION["Equirectangular"],',
  'PROJECTION["Mercator_2SP"],',
  #conic
  'PROJECTION["Albers_Conic_Equal_Area"],',
  'PROJECTION["Equidistant_Conic"],',  
  'PROJECTION["Lambert_Conformal_Conic_2SP"],',
  #transverse
  'PROJECTION["Transverse_Mercator"],',
  #oblique
  'PROJECTION["Hotine_Oblique_Mercator"],'
)
lat <- 40
lon <- -85
azi <- 335
rga <- azi
scf <- 0.9996
feast <- 0
fnorth <- 0
stp1 <- lat
stp2 <- lat+8
lorg <- lat-8
cmer <- lon

par1 <- paste0('PARAMETER["latitude_of_center",',lat,'],')
par2 <- paste0('PARAMETER["longitude_of_center",',lon,'],')
par3 <- paste0('PARAMETER["false_easting",',feast,'],')
par4 <- paste0('PARAMETER["false_northing",',fnorth,'],')
par5 <- paste0('PARAMETER["scale_factor",',scf,'],')
par6a <- paste0('PARAMETER["azimuth",',azi,'],')
par6b <- paste0('PARAMETER["rectified_grid_angle",',rga,'],')

par1a <- paste0('PARAMETER["standard_parallel_1",',stp1,'],')
par1b <- paste0('PARAMETER["standard_parallel_2",',stp2,'],')
par1c <- paste0('PARAMETER["latitude_of_origin",',lorg,'],')
par2a <- paste0('PARAMETER["central_meridian",',cmer,'],')



ifaz <- paste0(begining,
               datums[1],
               projections[1],
               par1,par2,par3,par4,units[1],axises
)
ifcyl <- paste0(begining,
                datums[1],
                projections[4],
                par1a,par2,par3,par4,units[1],axises
)
ifcon <- paste0(begining,
                datums[1],
                projections[4],
                par1,par2,par1a,par1b,par3,par4,units[1],axises
)
ifconf <- paste0(begining,
                 datums[1],
                 projections[4],
                 par1,par2a,par1a,par1b,par3,par4,units[1],axises
)
ifutm <- paste0(begining,
                datums[1],
                projections[10],
                par1,par2a,par5,par3,par4,units[1],axises
)

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
'PROJECTION["Hotine_Oblique_Mercator"],
    PARAMETER["latitude_of_center",45.3091666666667],
    PARAMETER["longitude_of_center",-86],
    PARAMETER["azimuth",337.25556],
    PARAMETER["rectified_grid_angle",337.25556],
    PARAMETER["scale_factor",0.9996],
    PARAMETER["false_easting",2546731.496],
    PARAMETER["false_northing",-4354009.816],
    UNIT["metre",1,
        AUTHORITY["EPSG","9001"]],
    AXIS["Easting",EAST],
    AXIS["Northing",NORTH]]'

CS[Cartesian,2],
AXIS["easting (X)",east,
     ORDER[1],
     LENGTHUNIT["metre",1]],
AXIS["northing (Y)",north,
     ORDER[2],
     LENGTHUNIT["metre",1]],
USAGE[
  SCOPE["State-wide spatial data management."],
  AREA["United States (USA) - Michigan."],
  BBOX[41.69,-90.42,48.32,-82.13]]'


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