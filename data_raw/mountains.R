library(climatools)
library(terra)
library(sf)

path <- 'C:/workspace2/dem/input'
input <- list.files(path)
mts <- read.csv('C:/workspace2/ClimateClassification/terrain/statemountains.csv')
mts <- st_read('C:/workspace2/climatools/data_raw/mountains.shp')
mts <- mts |> mutate(lat = st_coordinates(mts)[,2], lon = st_coordinates(mts)[,1])

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




for(k in 1:nrow(mts)){#k=9
thismts <- mts[k,]
cropdeg = 5
thisfile <- subset(df, west <= thismts$Longitude &
                     east >= thismts$Longitude &
                     north >= thismts$Latitude &
                     south <= thismts$Latitude)
lefttop <- subset(df, west <= thismts$Longitude-cropdeg/2 &
                     east >= thismts$Longitude-cropdeg/2 &
                     north >= thismts$Latitude+cropdeg/2 &
                     south <= thismts$Latitude+cropdeg/2)
righttop <- subset(df, west <= thismts$Longitude+cropdeg/2 &
                      east >= thismts$Longitude+cropdeg/2 &
                      north >= thismts$Latitude+cropdeg/2 &
                      south <= thismts$Latitude+cropdeg/2)
leftbottom <- subset(df, west <= thismts$Longitude-cropdeg/2 &
                     east >= thismts$Longitude-cropdeg/2 &
                     north >= thismts$Latitude-cropdeg/2 &
                     south <= thismts$Latitude-cropdeg/2)
rightbottom <- subset(df, west <= thismts$Longitude+cropdeg/2 &
                      east >= thismts$Longitude+cropdeg/2 &
                      north >= thismts$Latitude-cropdeg/2 &
                      south <= thismts$Latitude-cropdeg/2)
dem <- rast(paste0(path,'/',thisfile$filename))
demlist <- unique(c(paste0(path,'/',thisfile$filename),
paste0(path,'/',lefttop$filename),
paste0(path,'/',righttop$filename),
paste0(path,'/',leftbottom$filename),
paste0(path,'/',rightbottom$filename)))
for(i in 1:length(demlist)){
  dem0 <- rast(demlist[i])
  dem0 <- crop(dem0,ext(thismts$Longitude-cropdeg/2, thismts$Longitude+cropdeg/2, thismts$Latitude-cropdeg/2, thismts$Latitude+cropdeg/2))
  if(i==1){dem <- dem0}else{dem <- merge(dem,dem0)}
}
plot(dem)
rg <- minmax(dem)[2]-minmax(dem)[1]
scl <- rg*15+10000
demx <- reproject(dem=dem, rs=250, w=scl, h=scl, lat=thismts$Latitude, lon=thismts$Longitude)

plot(demx)
prebreaks = c(150,300,500,1000,2000,3000,4000,5000,6000,7000,8000)
maxrange <- max(minmax(demx)[2]-minmax(demx)[1],1000)
breaks <- prebreaks[prebreaks <=maxrange]
rng <- getrelief(demx, r1=1000, r2=maxrange, s=0.1, n=1, p='medium', breaks = breaks)
# plot(rng)


rng2 <- getrelief(demx, r1=1000, r2=maxrange, s=0.25, n=1, p='medium', breaks = breaks)
# plot(rng2)

steepslope = minmax(rng2)[2]
totalslope = minmax(rng)[2]
minmax(demx)[2] - minmax(demx)[1]
summit = minmax(demx)[2]

summitcoord = getmaxcoords(demx)
totalcoord = getmaxcoords(rng)
steepcoord = getmaxcoords(rng2)

mountains0 <- data.frame(state=mts$State[k], Name=mts$Name[k], summit=summit,totalslope=totalslope,steepslope=steepslope,summitcoord=summitcoord,totalcoord=totalcoord,steepcoord=steepcoord)
if(k==1){mountains=mountains0}else{mountains=rbind(mountains,mountains0)}

}
write.csv(mountains, 'C:/workspace2/dem/output/mountains.csv', row.names = F)


sl<- terrain(demx, v="slope", unit="radians")
as<- terrain(demx, v="aspect", unit="radians")
hsd <- terra::shade(sl,as)
plot(hsd)
writeRaster(hsd,'C:/workspace2/dem/output/orizaba/hsd.tif')
writeRaster(rng,'C:/workspace2/dem/output/orizaba/rng.tif')
writeRaster(demx,'C:/workspace2/dem/output/orizaba/dem.tif')



#----test----
prebreaks = c(150,300,500,1000,2000,3000,4000,5000,6000,7000,8000)
maxrange <- ceiling(max(minmax(demx)[2]-minmax(demx)[1],150))
breaks <- prebreaks[prebreaks <=maxrange]
rng <- getrelief(demx, r1=1000, r2=maxrange, s=0.1, n=1, p='medium', breaks = breaks)
plot(rng)
writeRaster(rng,paste0('C:/workspace2/dem/output/test/','rng','.tif'), overwrite=T)
bks <- c(1000, breaks*10, maxrange*10)
for(i in 1:length(bks)){
rng <- getrelief(demx, r1=bks[i], s=0.1, n=0, p='medium')
names(rng) <- paste0('r', bks[i])
plot(rng)
writeRaster(rng,paste0('C:/workspace2/dem/output/test/',paste0('r', bks[i]),'.tif'), overwrite=T)
}
path <- 'C:/workspace2/dem/output/test'
input <- list.files(path)
nms <- (stringr::str_split_fixed(input, '\\.', n=2)[,1])
for (i in 1:length(input)){#i=1
  x = rast(paste0(path,'/',input[i]))
    assign(paste0(nms[i]), x)
}
rngstack <- rast(mget(nms))
getsamp <- terra::spatSample(rngstack, 15)

xrng <- getsamp[11,'rng']
x2 <- getsamp[11,'r10000']
x1 <- getsamp[11,'r5000']
x2
x2/10000
x1
x1/5000
xrng
x2*(x1/5000 - 0.1)/(x1/5000-x2/10000)+x1*(0.1 - x2/10000)/(x1/5000-x2/10000)
