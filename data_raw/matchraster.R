library(terra)
library(climatools)
altlayer <- rast(system.file("extdata", "altlayer.tif", package="climatools") )
climatedata <- read.csv(system.file("extdata", "pts.p.csv", package="climatools"))
pts <- data.frame(x=climatedata$lon, y=climatedata$lat, z=climatedata$p07)
pts$z <- log10(pts$z + 1)
cropto <- c(-126,-105,35,50)
p07 <- toclimrast(pts, altlayer, cropto, covrange = 500, segx = 10, segy = 10); names(p07) <- 'p07'
plot(p07, col=rev(map.pal('bcyr')))
r1 <- 10^p07-1

r2 <- refitrast(p07mm, altlayer$elev, rotations=3)




plot(r1-r2)

r1max <- focal(r1, w=3, fun = 'max')
r1min <- focal(r1, w=3, fun = 'min')
r1med <- focal(r1, w=3, fun = 'median')
r2max <- focal(r2, w=3, fun = 'max')
r2min <- focal(r2, w=3, fun = 'min')
r2med <- focal(r2, w=3, fun = 'median')
r3 <- (r1med - r2med) + r2
r1top <- r1max-r1med; r1top <- ifel(r1top > 0,r1top,0)
r2top <- r2max-r2med; r2top <- ifel(r2top > 0,r2top,0)
r1bot <- r1min-r1med; r1bot <- ifel(r1bot < 0,r1bot,0)
r2bot <- r2min-r2med; r2bot <- ifel(r2bot < 0,r2bot,0)
plot(r2bot)

