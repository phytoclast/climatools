library(climatools)
library(terra)
library(sf)

path <- 'C:/workspace2/dem/input'
input <- list.files(path)
mts <- read.csv('C:/workspace2/ClimateClassification/terrain/statemountains.csv')
mts <- st_read('C:/workspace2/climatools/data_raw/mountains.shp')
mts <- readRDS('C:/workspace2/climatools/data_raw/mountains.RDS')
mts <- mts |> mutate(lat = st_coordinates(mts)[,2], lon = st_coordinates(mts)[,1])
rownames(mts) <- 1:nrow(mts)
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





k=49
thismts <- mts[k,]

cropdeg = 6
thisfile <- subset(df, west <= thismts$lon &
                     east >= thismts$lon &
                     north >= thismts$lat &
                     south <= thismts$lat)
lefttop <- subset(df, west <= thismts$lon-cropdeg/2 &
                    east >= thismts$lon-cropdeg/2 &
                    north >= thismts$lat+cropdeg/2 &
                    south <= thismts$lat+cropdeg/2)
righttop <- subset(df, west <= thismts$lon+cropdeg/2 &
                     east >= thismts$lon+cropdeg/2 &
                     north >= thismts$lat+cropdeg/2 &
                     south <= thismts$lat+cropdeg/2)
leftbottom <- subset(df, west <= thismts$lon-cropdeg/2 &
                       east >= thismts$lon-cropdeg/2 &
                       north >= thismts$lat-cropdeg/2 &
                       south <= thismts$lat-cropdeg/2)
rightbottom <- subset(df, west <= thismts$lon+cropdeg/2 &
                        east >= thismts$lon+cropdeg/2 &
                        north >= thismts$lat-cropdeg/2 &
                        south <= thismts$lat-cropdeg/2)
dem <- rast(paste0(path,'/',thisfile$filename))
demlist <- unique(c(paste0(path,'/',thisfile$filename),
                    paste0(path,'/',lefttop$filename),
                    paste0(path,'/',righttop$filename),
                    paste0(path,'/',leftbottom$filename),
                    paste0(path,'/',rightbottom$filename)))
demlist <- demlist[grepl('\\.tif$',demlist)]
for(i in 1:length(demlist)){#i=2
  dem0 <- rast(demlist[i])
  dem0 <- crop(dem0,ext(thismts$lon-cropdeg/2, thismts$lon+cropdeg/2, thismts$lat-cropdeg/2, thismts$lat+cropdeg/2))
  if(i==1){demy <- dem0}else{demy <- merge(demy,dem0)}
}
# plot(demy);plot(st_geometry(mts), add=T)

dem0 <- (crop(demy, ext(-153, -149, 62.5, 63.5)))
writeRaster(dem0, 'C:/workspace2/climatools/data_raw/denali.tif', overwrite=TRUE)
#reproject
# rg <- minmax(demy)[2]-minmax(demy)[1]
rg <- thismts$ht-minmax(demy)[1]

scl <- rg*2*10+10000
demx0 <- reproject(dem=demy, rs=250, w=100000, h=100000, lat=thismts$lat, lon=thismts$lon, method='bilinear')
demx <- climatools::RestoreMaxMin(demy,demx0,mts,'ht')
plot(dem0)
prebreaks = c(50, 100,150,200,300,500,1000,2000,3000,4000,5000,6000,7000,8000)
maxrange <- max(rg2,1000)
breaks <- prebreaks[prebreaks <=maxrange]
rng <- getrelief(demx, r1=1000, r2=maxrange, s=0.1, n=1, p='medium', breaks = breaks)
plot(rng)
hsd <- hillshade(demx)
plot(hsd)
rng2 <- getrelief(demx, r1=1000, r2=maxrange, s=0.25, n=1, p='medium', breaks = breaks)
plot(rng2)
names(demx) <- 'elevation'
names(hsd) <- 'hillshade'
names(rng) <- 'broadrelief'
names(rng2) <- 'steeprelief'
rastcomp <- c(demx,hsd,rng,rng2)
writeRaster(rastcomp, 'C:/workspace2/dem/output/rastcomp.tif', overwrite=TRUE)






#maximum circle
mtcrop <- thismts |> vect() |> terra::project(demx) |> crop(ext(demx))
mtbuff <- mtcrop |> buffer(pmax(width=rg*10.1, 2020))
demx <- demx |> crop(mtbuff)
mtzone <- mtbuff |> rasterize(demx) |> crop(ext(mtbuff))
# plot(demx);plot(mtbuff, add=T)

#get range within big circle
crange <- minmax(demx*mtzone)
rg2 <- thismts$ht-crange[1]
cutoff <- (thismts$ht-crange[1])*0.5+crange[1]
mtbuff2 <- mtcrop |> buffer(pmax(width=rg2*5.05, 1010))
mtzone2 <- mtbuff2 |> rasterize(demx)
mtbuff2.5 <- mtcrop |> buffer(pmax(width=rg2*1.25, 1010))
mtzone2.5 <- mtbuff2.5 |> rasterize(demx)
#shave off tops of adjacent mountains
demxcut <- ifel(demx > cutoff  & is.na(mtzone2.5), cutoff, demx)
# plot(demxcut)
#get 10% range within mountain's zone
prebreaks = c(50, 100,150,200,300,500,1000,2000,3000,4000,5000,6000,7000,8000)
maxrange <- max(rg2,1000)
breaks <- prebreaks[prebreaks <=maxrange]
rng <- getrelief(demxcut, r1=1000, r2=maxrange, s=0.1, n=1, p='medium', breaks = breaks)
# plot(rng)


#restrict range to near mountain
rng <- mtzone2*rng
# plot(rng)
newrg <- minmax(rng)[2]
mtbuff3 <- mtcrop |> buffer(pmax(width=newrg*5.05, 1010))
mtzone3 <- mtbuff3 |> rasterize(demx)
rng <- mtzone3*rng

#25% relief
rng2 <- getrelief(demxcut, r1=1000, r2=maxrange, s=0.25, n=1, p='medium', breaks = breaks)
# plot(rng2)
rng2 <- mtzone3*rng2
# plot(rng2)


broadrelief = minmax(rng)[2]
steeprelief = minmax(rng2)[2]
# summit = minmax(demx)[2]

# summitcoord = getmaxcoords(demx)
broadcoord = getmaxcoords(rng)
steepcoord = getmaxcoords(rng2)








temp <- rast('C:/workspace2/processclimategrids/output/Tw.tif')
denali0 <- rast('data_raw/denali.tif')
Tw <- crop(temp, denali0)
plot(Tw)
writeRaster(Tw, 'C:/workspace2/climatools/data_raw/Tw.tif', overwrite=T)


Tw0 <- rast('data_raw/Tw.tif')
Tw <- list(data = as.matrix(Tw0,wide=TRUE), crs = crs(Tw0), ext = t(as.matrix(ext(Tw0))), names = "Temperature.WarmMonth.Worldclim2")
usethis::use_data(Tw, overwrite = T)

denali0 <- rast('data_raw/denali.tif')
denali <- list(data = as.matrix(denali0,wide=TRUE),
               crs = crs(denali0),
               ext = t(as.matrix(ext(denali0))),
               names = names(denali0))
usethis::use_data(denali, overwrite = T)

peaks <- read.csv('C:/workspace2/climatools/data_raw/worldpeaks4.csv')
peaks <- peaks[,c("Continent","state","Name","summit","broadrelief","steeprelief","lat","lon")]
colnames(peaks) <- c("continent","state","name","summit","relief.broad","relief.steep","lat","lon")
peaks <- peaks |> mutate(x=lon,y=lat)
peaks <- st_as_sf(peaks, coords = c("x","y"), crs = 4326)
usethis::use_data(peaks, overwrite = T)

fromraster <- function(x){
  y <- list(data = as.matrix(x,wide=TRUE),
            crs = crs(x),
            ext = t(as.matrix(ext(x))),
            names = names(x))
  return(y)}

toraster <- function(x){
  y <- rast(x$data, crs = x$crs, ext = ext(x$ext))
  names(y) <- x$names
  return(y)}


data("denali")
denali <- denali
denali <- toraster(denali)
minmax(denali)
#reproject to a coarser resolution
denali0 <- reproject(denali, lat = 63, lon = -151, rs = 1000,  h=100000, w=200000)
#note values less extreme than original
minmax(denali0)
denali1 <- RestoreMaxMin(denali, denali0)
#note max/min values are restored
minmax(denali1)

#supplemental point elevations
data("peaks")
peaks <- peaks
#reproject to a finer resolution
denali0 <- reproject(denali, lat = 63, lon = -151, rs = 100,  h=100000, w=200000)
peaks <- st_transform(peaks, crs(denali0))
plot(denali0); plot(st_geometry(peaks), add=TRUE)
#Note that highest peak is known to be higher than raster indicates even if preserving original raster values when projected at higher resolution.
peaks[peaks$name %in% 'Denali',]$summit
minmax(denali0)
denali1 <- RestoreMaxMin(x=denali, y=denali0, s=peaks, e='summit')
#New high point is equal to supplemental point data.
minmax(denali1)

hsd <- hillshade(denali0)
plot(hsd)

#load package data for elevation
data("denali")
dem0 <- denali
dem0 <- toraster(dem0)
dem <- reproject(dem0, rs=250)
plot(dem)

#load package data for July temperature
data("Tw")
temperature0 <- Tw
temperature0 <- toraster(temperature0)

#Original temperature is coarse resolution.
plot(temperature0)

#This resamples to match resolution of elevation, but it appears blurry.
temperature <- project(temperature0, dem)
plot(temperature)

#prepare lower resolution elevation raster matching temperature resolution
dem0 <- project(dem0, temperature0)
plot(dem0)

#This estimates a local relationship of temperature with elevation.
temperature <- enhanceRast(temperature0, dem0,  dem)
plot(temperature)


data(Tw)
Tw0 <- Tw
Tw <- rast(Tw0$dem, crs = Tw0$crs, ext = ext(Tw0$ext))
names(denali) <- denali0$names; rm(denali0)
plot(denali)



hsd <- hillshade(denali)
plot(hsd)
dem <- reproject(dem=denali, rs=250, w=100000, h=50000, lat= 63, lon= -151, method='bilinear')
peaks2 <- st_transform(peaks, crs=crs(dem))
plot(dem); plot(st_geometry(peaks2), add=TRUE)

temperature <- data(Tw)
plot(temperature)
demcoarse <- terra::project(denali, Tw, method='bilinear')
enhanceRast()


#load package data for elevation
data("denali")
dem <- denali
dem <- toraster(dem)
plot(dem)
#Reproject to get metric distance units.
dem <- reproject(dem, rs=250)

#Established desired vertical increment of focal range analysis.
prebreaks = c(50, 100,150,200,300,500,1000,2000,3000,4000,5000,6000,7000,8000)
maxrange <- max(minmax(dem)[2],1000)#Establish maximum relative to the highest elevation in DEM. This is to reduce the need for radii wider than required.
breaks <- prebreaks[prebreaks <=maxrange]
#Relief at a 10% slope. Using "low" precision for faster performance
rng10 <- getrelief(dem, r1=1000, r2=maxrange*10, s=0.1, n=1, p='low', breaks = breaks)
plot(rng10)
#Relief at a 25% slope. Using "medium" precision for better quality.
rng25 <- getrelief(dem, r1=1000, r2=maxrange*10, s=0.25, n=1, p='medium', breaks = breaks)
plot(rng25)



