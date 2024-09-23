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


thislon <- -100
thislat <- 45


cropdeg = 20
thisfile <- subset(df, west <= thislon &
                     east >= thislon &
                     north >= thislat &
                     south <= thislat)
lefttop <- subset(df, west <= thislon-cropdeg/2 &
                    east >= thislon-cropdeg/2 &
                    north >= thislat+cropdeg/2 &
                    south <= thislat+cropdeg/2)
righttop <- subset(df, west <= thislon+cropdeg/2 &
                     east >= thislon+cropdeg/2 &
                     north >= thislat+cropdeg/2 &
                     south <= thislat+cropdeg/2)
leftbottom <- subset(df, west <= thislon-cropdeg/2 &
                       east >= thislon-cropdeg/2 &
                       north >= thislat-cropdeg/2 &
                       south <= thislat-cropdeg/2)
rightbottom <- subset(df, west <= thislon+cropdeg/2 &
                        east >= thislon+cropdeg/2 &
                        north >= thislat-cropdeg/2 &
                        south <= thislat-cropdeg/2)
dem <- rast(paste0(path,'/',thisfile$filename))
demlist <- unique(c(paste0(path,'/',thisfile$filename),
                    paste0(path,'/',lefttop$filename),
                    paste0(path,'/',righttop$filename),
                    paste0(path,'/',leftbottom$filename),
                    paste0(path,'/',rightbottom$filename)))
demlist <- demlist[grepl('\\.tif$',demlist)]
for(i in 1:length(demlist)){#i=2
  dem0 <- rast(demlist[i])
  dem0 <- crop(dem0,ext(thislon-cropdeg/2, thislon+cropdeg/2, thislat-cropdeg/2, thislat+cropdeg/2))
  if(i==1){demy <- dem0}else{demy <- merge(demy,dem0)}
}
 plot(dem);plot(st_geometry(mts), add=T)

 dem <- rast(paste0(path,'/',df$filename[28]))
 prj = setProjection(prj = 'equalarea.conic', 
                     lat=mean(ext(dem)[3:4])-7, 
                     lon=mean(ext(dem)[1:2]),
                     lat2=mean(ext(dem)[3:4])+7)
 demx <- reproject(dem, prj=prj, rs = 500)
 plot(demx);plot(st_geometry(st_transform(mts,crs(demx))), add=T)
 
 demx <- RestoreMaxMin(dem,demx, mts,'ht')
 prebreaks = c(300,500,1000,1500,2000,3000,4000,5000,6000,7000,8000)
 maxrange <- max(minmax(demx)[2]-minmax(demx)[1],10000)
 breaks <- prebreaks[prebreaks <=maxrange]
 rng <- getrelief(demx, r1=3000, r2=maxrange, s=0.1, n=1, p='medium', breaks = breaks)
 plot(rng)
 
 
 (49-29)*
(45-29)/(49-26)