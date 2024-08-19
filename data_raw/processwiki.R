library(climatools)
proms <- read.delim('C:/workspace2/dem/prom1.txt')
proms2 <- read.delim('C:/workspace2/dem/usprompeaks.txt')
proms3 <- read.delim('C:/workspace2/dem/statepeaks.txt')
proms4 <- read.delim('C:/workspace2/dem/otherpeaks.txt')
proms5 <- read.delim('C:/workspace2/dem/isolatedpeaks.txt')
colnames(proms) <- c("No","Peak","Range","Location",
                     "Coords","Prominence","ht","Col",
                     "parent", "Prom.parent")
colnames(proms2) <- c("Rank","Peak","State","Range","ht","Prominence","Isolation","Coords" )
colnames(proms3) <- c("State", "Peak","ht","Coords" )
colnames(proms4) <- c("Peak", "State","ht","Coords" )
colnames(proms5) <- c("Peak", "Region","Range","ht","Prominence","Isolation", "Coords" )

proms <- proms |> mutate(lat0 = stringr::str_split_fixed(Coords, 'N|S',2)[,1],
                         lon0 = stringr::str_split_fixed(Coords, 'N|S',2)[,2],
                         lat1 = stringr::str_split_fixed(lat0, '°|′|″',4)[,1],
                         lat2 = stringr::str_split_fixed(lat0, '°|′|″',4)[,2],
                         lat3 = stringr::str_split_fixed(lat0, '°|′|″',4)[,3],
                         lon1 = stringr::str_split_fixed(lon0, '°|′|″|E|W',4)[,1],
                         lon2 = stringr::str_split_fixed(lon0, '°|′|″|E|W',4)[,2],
                         lon3 = stringr::str_split_fixed(lon0, '°|′|″|E|W',4)[,3],
                         lon1 = stringr::str_extract(lon1,'[:digit:]+'),
                         lat1 = as.numeric(trimws(lat1)),
                         lat2 = as.numeric(lat2),
                         lat3 = as.numeric(lat3),lat3 = ifelse(is.na(lat3),0,lat3),
                         lon1 = as.numeric(trimws(lon1)),
                         lon2 = as.numeric(lon2),
                         lon3 = as.numeric(lon3),lon3 = ifelse(is.na(lon3),0,lon3),
                         
                         lat = lat1+lat2/60+lat3/3600,lon = lon1+lon2/60+lon3/3600,
                         lat = ifelse(grepl('S',Coords),lat*-1,lat),
                         lon = ifelse(grepl('W',Coords),lon*-1,lon))

proms2 <- proms2 |> mutate(lat0 = stringr::str_split_fixed(Coords, 'N|S',2)[,1],
                         lon0 = stringr::str_split_fixed(Coords, 'N|S',2)[,2],

                         lat = stringr::str_extract(lat0,'[:digit:]+\\.[:digit:]*'),
                         lon = stringr::str_extract(lon0,'[:digit:]+\\.[:digit:]*'),
                         lat = as.numeric(lat),lon = as.numeric(lon),
                         lat = ifelse(grepl('S',Coords),lat*-1,lat),
                         lon = ifelse(grepl('W',Coords),lon*-1,lon),
                         Peak=stringr::str_split_fixed(Peak, '\\[',2)[,1],
                         ht0 = stringr::str_split_fixed(ht, 'ft',3)[,2],
                         ht01 = stringr::str_split_fixed(ht, 'm',3)[,1],
                         ht0 = ifelse(nchar(ht0)<=1, ht01,ht0),
                         ht1 = stringr::str_extract(ht0,'[:digit:]+\\.*[:digit:]*'),
                         ht = as.numeric(ht1))

proms3 <- proms3 |> mutate(lat0 = stringr::str_split_fixed(Coords, 'N|S',2)[,1],
                           lon0 = stringr::str_split_fixed(Coords, 'N|S',2)[,2],
                           
                           lat = stringr::str_extract(lat0,'[:digit:]+\\.[:digit:]*'),
                           lon = stringr::str_extract(lon0,'[:digit:]+\\.[:digit:]*'),
                           lat = as.numeric(lat),lon = as.numeric(lon),
                           lat = ifelse(grepl('S',Coords),lat*-1,lat),
                           lon = ifelse(grepl('W',Coords),lon*-1,lon),
                           Peak=stringr::str_split_fixed(Peak, '\\[',2)[,1],
                           State=stringr::str_split_fixed(State, '\\[',2)[,1],
                           ht1 = stringr::str_extract(ht,'[:digit:]+\\.*[:digit:]*'),
                           ht = as.numeric(ht1))

proms4 <- proms4 |> mutate(lat0 = stringr::str_split_fixed(Coords, 'N|S',2)[,1],
                           lon0 = stringr::str_split_fixed(Coords, 'N|S',2)[,2],
                           lat1 = stringr::str_split_fixed(lat0, '°|′|\'|\'\'|″',4)[,1],
                           lat2 = stringr::str_split_fixed(lat0, '°|′|\'|\'\'|″',4)[,2],
                           lat3 = stringr::str_split_fixed(lat0, '°|′|\'|\'\'|″',4)[,3],
                           lon1 = stringr::str_split_fixed(lon0, '°|′|\'|\'\'|″|E|W',4)[,1],
                           lon2 = stringr::str_split_fixed(lon0, '°|′|\'|\'\'|″|E|W',4)[,2],
                           lon3 = stringr::str_split_fixed(lon0, '°|′|\'|\'\'|″|E|W',4)[,3],
                           lon1 = stringr::str_extract(lon1,'[:digit:]+'),
                           
                           lat1 = as.numeric(trimws(lat1)),
                           lat2 = as.numeric(lat2),
                           lat3 = as.numeric(lat3),lat3 = ifelse(is.na(lat3),0,lat3),
                           lon1 = as.numeric(trimws(lon1)),
                           lon2 = as.numeric(lon2),
                           lon3 = as.numeric(lon3),lon3 = ifelse(is.na(lon3),0,lon3),
                           
                           lat = lat1+lat2/60+lat3/3600,lon = lon1+lon2/60+lon3/3600,
                           lat = ifelse(grepl('S',Coords),lat*-1,lat),
                           lon = ifelse(grepl('W',Coords),lon*-1,lon),
                           Peak=stringr::str_split_fixed(Peak, '\\[',2)[,1],
                           State=stringr::str_split_fixed(State, '\\[',2)[,1],
                           ht1 = stringr::str_extract(ht,'[:digit:]+\\.*[:digit:]*'),
                           ht = as.numeric(ht1))

proms5 <- proms5 |> mutate(lat0 = stringr::str_split_fixed(Coords, 'N|S',2)[,1],
                           lon0 = stringr::str_split_fixed(Coords, 'N|S',2)[,2],
                           
                           lat = stringr::str_extract(lat0,'[:digit:]+\\.[:digit:]*'),
                           lon = stringr::str_extract(lon0,'[:digit:]+\\.[:digit:]*'),
                           lat = as.numeric(lat),lon = as.numeric(lon),
                           lat = ifelse(grepl('S',Coords),lat*-1,lat),
                           lon = ifelse(grepl('W',Coords),lon*-1,lon),
                           Peak=stringr::str_split_fixed(Peak, '\\[',2)[,1],
                           ht0 = stringr::str_split_fixed(ht, 'ft',3)[,2],
                           ht01 = stringr::str_split_fixed(ht, 'm',3)[,1],
                           ht0 = ifelse(nchar(ht0)<=1, ht01,ht0),
                           ht1 = stringr::str_extract(ht0,'[:digit:]+\\.*[:digit:]*'),
                           ht = as.numeric(ht1))

proms <- subset(proms, select=c(Peak, Range, lat, lon, ht))
proms2 <- subset(proms2, select=c(Peak, Range, lat, lon, ht))
proms3 <- subset(proms3, select=c(Peak, State, lat, lon, ht))
proms4 <- subset(proms4, select=c(Peak, State, lat, lon, ht))
proms5 <- subset(proms5, select=c(Peak, Range, lat, lon, ht))
promx <- bind_rows(proms, proms2, proms3, proms4,proms5)
library(sf)
states <- st_read('C:/a/geo/world/level4/level4.shp')
promsf <- st_as_sf(promx, coords = c("lon", "lat"), crs = 4326)
st_crs(states) <- st_crs(4326); states <- st_make_valid(states);# plot(st_geometry(states))
promsf <- st_join(promsf, states)
promsf <- promsf |> mutate(State=ifelse(is.na(State), LEVEL_4_NA, State))
promsf <- promsf |> subset(select=c(Peak,Range,State,CONTINENT, ht))
st_write(promsf, dsn = "output/corrected", layer = "mountains.shp", driver = 'ESRI Shapefile', append = F)
#rasters
library(terra)
path <- 'C:/workspace2/dem/input'
input <- list.files(path)
i=5
dem0 <- rast(paste0(path,'/',input[i]))
plot(dem0)
plot(st_geometry(promsf), add=T)
90/10000000*1000
demax <- focalmax(dem0, r=0.009)
fake <- rast(xmin=ext(dem0)[1], xmax=ext(dem0)[2],
             ymin=ext(dem0)[3], ymax=ext(dem0)[4], crs=crs(dem0), res=res(dem0)*4)
dotrast <- rasterize(promsf,fake, field="ht")
dotrast <- project(dotrast, demax, method='near')
dotrast0 <- focalmax(dotrast, r=0.009)

newrast <- ifel(!is.na(dotrast0) & dotrast > dem0 & demax == dem0, dotrast0, dem0,
                filename=paste0("output/corrected/", input[i]), overwrite=T)


