library(climatools)
library(terra)
library(sf)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
path <- 'C:/scripts/dem/data'
input <- list.files(path)
mts <- climatools::peaks
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

llats <- c(0:10)/10*(80-10)+10
llons <- c(0:10)/10*(-45--179)+-180
ll <- merge(llats, llons); colnames(ll) <- c('lat','lon')
(80-10)/20
(-30--180)/20
for(k in 1:nrow(ll)){
#k=220
thismts <- ll[k,]

croplat = (llats[2]-llats[1])+2
croplon = (llons[2]-llons[1])+2/cos((thismts$lat)/360*2*3.141592)
croplat2 = (llats[2]-llats[1])+1
croplon2 = (llons[2]-llons[1])+1/cos((thismts$lat)/360*2*3.141592)

thisfile <- subset(df, west <= thismts$lon &
                     east >= thismts$lon &
                     north >= thismts$lat &
                     south <= thismts$lat)
lefttop <- subset(df, west <= thismts$lon-croplon/2 &
                    east >= thismts$lon-croplon/2 &
                    north >= thismts$lat+croplat/2 &
                    south <= thismts$lat+croplat/2)
righttop <- subset(df, west <= thismts$lon+croplon/2 &
                     east >= thismts$lon+croplon/2 &
                     north >= thismts$lat+croplat/2 &
                     south <= thismts$lat+croplat/2)
leftbottom <- subset(df, west <= thismts$lon-croplon/2 &
                       east >= thismts$lon-croplon/2 &
                       north >= thismts$lat-croplat/2 &
                       south <= thismts$lat-croplat/2)
rightbottom <- subset(df, west <= thismts$lon+croplon/2 &
                        east >= thismts$lon+croplon/2 &
                        north >= thismts$lat-croplat/2 &
                        south <= thismts$lat-croplat/2)
if(nrow(thisfile) > 0){
dem <- rast(paste0(path,'/',thisfile$filename))
demlist <- unique(c(paste0(path,'/',thisfile$filename),
                    paste0(path,'/',lefttop$filename),
                    paste0(path,'/',righttop$filename),
                    paste0(path,'/',leftbottom$filename),
                    paste0(path,'/',rightbottom$filename)))
demlist <- demlist[grepl('\\.tif$',demlist)]
for(i in 1:length(demlist)){#i=2
  dem0 <- rast(demlist[i])
  dem0 <- crop(dem0,ext(thismts$lon-croplon/2, thismts$lon+croplon/2, thismts$lat-croplat/2, thismts$lat+croplat/2))
  if(i==1){demy <- dem0}else{demy <- merge(demy,dem0)}
}
if(max(minmax(demy))-min(minmax(demy)) > 0){
plot(demy);plot(st_geometry(mts), add=T)
demx0 <- reproject(dem=demy, rs=500,  lat=thismts$lat, lon=thismts$lon, method='bilinear')
demx1 <- climatools::RestoreMaxMin(demy,demx0,mts,'summit')
plot(demx1)
rvals <- c(200,300,500,700,1000,1500,2000,3000,4000,5000,6000,7000)
rng <- getrelief(demx1,  n=1,s=0.1,breaks = rvals, p='medium')
st_crs(demy)
rng2 <- reproject(rng, prj=crs(demy))
rng2 <- climatools::RestoreMaxMin(rng,rng2)
plot(rng2)
rngcrop <- crop(rng2,ext(thismts$lon-croplon2/2, thismts$lon+croplon2/2, thismts$lat-croplat2/2, thismts$lat+croplat2/2))
plot(rngcrop)
writeRaster(rngcrop, paste0('C:/scripts/dem/output/rng',k,'.tif'), overwrite=TRUE)
}}
}

#########
library(climatools)
library(terra)
library(sf)


path <- 'C:/scripts/dem/output'
input <- list.files(path)
input <- input[grepl('\\.tif$', input)]

for(i in 1:length(input)){
rng0 <- rast(paste0(path,'/',input[i]))
assign(input[i], rng0)

} 

rng <- sprc(mget(input))
rng <- merge(rng)
plot(rng)
names(rng) <- 'rng'
writeRaster(rng, 'C:/scripts/dem/rng.tif', overwrite=TRUE)