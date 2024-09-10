library(climatools)
library(terra)
library(sf)

path <- 'C:/workspace2/dem/input'
input <- list.files(path)
mts <- read.csv('C:/workspace2/ClimateClassification/terrain/statemountains.csv')
mts <- st_read('C:/workspace2/climatools/data_raw/mountains.shp')
mts <- readRDS('C:/workspace2/climatools/data_raw/mountains.RDS')
mts <- mts |> mutate(lat = st_coordinates(mts)[,2], lon = st_coordinates(mts)[,1])

rownames(mts) <- NULL
getmaxcoords <- function(x){
coord <- as.data.frame(xyFromCell(x, (where.max(x)[2])))
coord <- st_as_sf(coord, coords = c("x","y"), crs = crs(x))
coord <- st_transform(coord, crs = 4326)
coord <- st_coordinates(coord)
return(coord)}

# mts2 <- mountains |> subset(summit < broadrelief)
# mts2 <- read.csv('C:/workspace2/dem/output/worldpeaks.csv')
# mts3 <- subset(mts, Peak %in% mts2$Name )
for(i in 1:length(input)){#i=1

  dem <- rast(paste0(path,'/',input[i]))
  df0 <- data.frame(i=i, filename=input[i], west=ext(dem)[1],east=ext(dem)[2],north=ext(dem)[4],south=ext(dem)[3])
  rownames(df0)=i
if(i==1){df=df0}else{df=rbind(df,df0)}
}



for(k in 1:nrow(mts)){#k=108
thismts <- mts[k,]
if(thismts$lat > -60){
cropdeg = 5
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



#reproject
# rg <- minmax(demy)[2]-minmax(demy)[1]
rg <- thismts$ht-minmax(demy)[1]
if(rg > 0){
scl <- rg*2*10+10000
demx0 <- reproject(dem=demy, rs=250, w=scl, h=scl, lat=thismts$lat, lon=thismts$lon, method='bilinear')
demx <- climatools::RestoreMaxMin(demy,demx0,mts,'ht')

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
mtbuff2.5 <- mtcrop |> buffer(pmax(width=rg2*3.05, 1010))
mtzone2.5 <- mtbuff2.5 |> rasterize(demx)
#shave off tops of adjacent mountains
demxcut <- ifel(demx > cutoff  & is.na(mtzone2.5), cutoff, demx)
# plot(demxcut)
#get 10% range within mountain's zone
prebreaks = c(50, 100,150,200,300,500,1000,1500,2000,3000,4000,5000,6000,7000,8000)
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


steeprelief = minmax(rng2)[2]
broadrelief = minmax(rng)[2]

# summit = minmax(demx)[2]

# summitcoord = getmaxcoords(demx)
broadcoord = getmaxcoords(rng)
steepcoord = getmaxcoords(rng2)

mountains0 <- data.frame(Continent=mts$CONTINENT[k], state=mts$State[k], Name=mts$Peak[k], summit=mts$ht[k],broadrelief=round(broadrelief,0),steeprelief=round(steeprelief,0),lat=mts$lat[k],lon=mts$lon[k],broadcoord=broadcoord,steepcoord=steepcoord)
if(k==1){mountains=mountains0}else{mountains=rbind(mountains,mountains0)}

}}}
write.csv(mountains, 'C:/workspace2/climatools/data_raw/worldpeaks4.csv', row.names = F)


