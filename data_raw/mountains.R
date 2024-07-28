library(climatools)
library(terra)
library(sf)

path <- 'C:/workspace2/dem/input'
input <- list.files(path)
mts <- read.csv('C:/workspace2/ClimateClassification/terrain/statemountains.csv')

getmaxcoords <- function(x){
coord <- as.data.frame(xyFromCell(x, (where.max(x)[2])))
coord <- st_as_sf(coord, coords = c("x","y"), crs = crs(x))
coord <- st_transform(coord, crs = 4326)
coord <- st_coordinates(coord)
return(coord)}


for(i in 1:length(input)){#i=1
  
  dem <- rast(paste0(path,'/',input[i]))
  df0 <- data.frame(i=i, filename=input[i], west=ext(dem)[1],east=ext(dem)[2],north=ext(dem)[4],south=ext(dem)[3])
  rownames(df0)=i
if(i==1){df=df0}else{df=rbind(df,df0)}
  
  }

thismts <- mts[59,]
cropdeg = 5
thisfile <- subset(df, west <= thismts$Longitude &
                     east >= thismts$Longitude &
                     north >= thismts$Latitude &
                     south <= thismts$Latitude)
leftfile <- subset(df, west <= thismts$Longitude-cropdeg/2 &
                     east >= thismts$Longitude-cropdeg/2 &
                     north >= thismts$Latitude &
                     south <= thismts$Latitude)
rightfile <- subset(df, west <= thismts$Longitude+cropdeg/2 &
                      east >= thismts$Longitude+cropdeg/2 &
                      north >= thismts$Latitude &
                      south <= thismts$Latitude)
bottomfile <- subset(df, west <= thismts$Longitude &
                     east >= thismts$Longitude &
                     north >= thismts$Latitude-cropdeg/2 &
                     south <= thismts$Latitude-cropdeg/2)
topfile <- subset(df, west <= thismts$Longitude &
                      east >= thismts$Longitude &
                      north >= thismts$Latitude+cropdeg/2 &
                      south <= thismts$Latitude+cropdeg/2)
dem <- rast(paste0(path,'/',thisfile$filename))
demlist <- unique(c(paste0(path,'/',thisfile$filename),
paste0(path,'/',leftfile$filename),
paste0(path,'/',rightfile$filename),
paste0(path,'/',bottomfile$filename),
paste0(path,'/',topfile$filename)))
for(i in 1:length(demlist)){
  dem0 <- rast(demlist[i])
  dem0 <- crop(dem0,ext(thismts$Longitude-cropdeg/2, thismts$Longitude+cropdeg/2, thismts$Latitude-cropdeg/2, thismts$Latitude+cropdeg/2))
  if(i==1){dem <- dem0}else{dem <- merge(dem,dem0)}
}
plot(dem)

demx <- reproject(dem=dem, rs=250, w=100000, h=100000, lat=thismts$Latitude, lon=thismts$Longitude)

plot(demx)

rng <- getrelief(demx, r1=1000, r2=90000, s=0.1, n=3, p='medium', breaks = c(150,300,1000,2500,5000,8000))
plot(rng);minmax(rng)[2]


rng2 <- getrelief(demx, r1=1000, r2=90000, s=0.25, n=3, p='medium', breaks = c(150,300,1000,2500,5000,7000))
plot(rng2);minmax(rng2)[2]
minmax(rng)[2]
minmax(rng2)[2]
minmax(demx)

getmaxcoords(rng)
getmaxcoords(rng2)
getmaxcoords(demx)


sl<- terrain(demx, v="slope", unit="radians")
as<- terrain(demx, v="aspect", unit="radians")
hsd <- terra::shade(sl,as)
plot(hsd)
writeRaster(hsd,'C:/workspace2/dem/output/orizaba/hsd.tif')
writeRaster(rng,'C:/workspace2/dem/output/orizaba/rng.tif')
writeRaster(demx,'C:/workspace2/dem/output/orizaba/dem.tif')









dem <- rast('C:/a/geo/dem/SRTM250m/usa/w001001.adf')
demx <- rast("C:/a/geo/dem/alt_30s_esri/alt/w001001.adf")
dem1 <-  reproject(dem=demx, w=500000, h=500000, lat=63, lon=-148)
plot(dem1)

dem2 <-  reproject(dem=dem, rs = 250, w=800000, h=750000, lat=38, lon=-108)
plot(dem2)


rng <- getrelief(dem1, r1=1000, r2=100000, s=0.1, n=9)
plot(rng)

rng2 <- (rng> 1000)+(rng> 2000)+(rng> 3000)+(rng> 4000)+(rng> 5000)+(rng> 6000)
plot(rng2)


dem2 <- rast('C:/workspace2/dem/50N150W_20101117_gmted_mea075.tif')
dem1 <- rast('C:/workspace2/dem/50N180W_20101117_gmted_mea075.tif')
dem <- merge(dem1,dem2)
plot(dem)

demx <-  reproject(dem=dem, rs = 250, w=800000, h=800000, lat=63, lon=-145)
plot(demx)
rng <- getrelief(demx, r1=1000, r2=100000, s=0.1, n=9)
plot(rng)
sl<- terrain(dem, v="slope", unit="radians")
as<- terrain(dem, v="aspect", unit="radians")
hsd <- terra::shade(sl,as)
plot(hsd)

writeRaster(dem, 'C:/workspace2/dem/output/dem.tif', overwrite=T)
writeRaster(hsd, 'C:/workspace2/dem/output/hsd.tif', overwrite=T)
writeRaster(rng, 'C:/workspace2/dem/output/rng.tif', overwrite=T)

demx <-  reproject(dem=dem, rs = 250)
rng <- getrelief(demx, r1=1000, r2=100000, s=0.1, n=9)
writeRaster(rng, 'C:/workspace2/dem/output/rng.tif', overwrite=T)


path <- 'C:/workspace2/dem/input'
input <- list.files(path)

# dem1 <- rast(paste0(path,'/',input[7]))
dem <- rast(paste0(path,'/',input[2]))
# dem <- merge(dem1,dem)
plot(dem)
demx <-  reproject(dem=dem, rs=250, w=500000, h=500000, lat=19.5, lon=-98.5)
plot(demx)
rng <- getrelief(demx, r1=1000, r2=100000, s=0.1, n=9)
plot(rng)
sl<- terrain(demx, v="slope", unit="radians")
as<- terrain(demx, v="aspect", unit="radians")
hsd <- terra::shade(sl,as)
plot(hsd)

writeRaster(hsd,'C:/workspace2/dem/output/orizaba/hsd.tif')
writeRaster(rng,'C:/workspace2/dem/output/orizaba/rng.tif')
writeRaster(demx,'C:/workspace2/dem/output/orizaba/dem.tif')




dem <- rast("C:/a/geo/dem/alt_30s_esri/alt/w001001.adf")
demx <-  reproject(dem=dem, w=100000, h=100000, lat=63, lon=-151)
plot(demx)


path <- 'C:/workspace2/dem/input'
input <- list.files(path)

dem <- rast(paste0(path,'/',input[11]))
# dem <- merge(dem1,dem)
plot(dem)
demy <-  reproject(dem=dem, w=100000, h=100000, lat=63, lon=-151)
plot(demy)


demz <- climatools::AmplifiedReduction(demy, fact=4)
plot(demz)
dema <-  reproject(dem=demz, w=100000, h=100000, lat=63, lon=-151, rs=250)
plot(dema)

x <- (1:10)/10-0.5
y <- x
df <- merge(x,y)
df$z <- 1000/((df$x-0.25)^2+(df$y-0.40)^2+1)^2
set.seed(1)
df$z <- df$z +rnorm(10*10, sd = 100)
df$x <- df$x*10 -85
df$y <- df$y*10 +20
dem <- rast(cbind(x=df$x,y=df$y,z=df$z), type="xyz", crs='EPSG:4326')
plot(dem)

dem <- reproject(dem=dem, rs = 250)
plot(dem)