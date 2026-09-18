library(climatools)
library(ranger)
dem <- climatools::denali |> climatools::toraster(); names(dem) <- 'elev'
t0 <- climatools::Tw |>  toraster()

dem0 <- project(dem, t0)
plot(dem0)
#This estimates a local relationship of temperature with elevation.
t1 <- enhanceRast(t0, dem0,  dem)
plot(t1)

library(sf)
x <- thisgrid#; t1<-logp(t1)
cropto <- c(-90, -75, 40, 55)
cropfrom <- c(-82,-80, 45.5,49)

modelfit <- function(x, elev=NULL, cropto, cropfrom=NULL, rotations=0){
  t0 <- crop(x, cropto); names(t0) <- 't0'
  xy0 <- makexyrast(t0,rotations)
  if(!is.null(cropfrom)){
    t00 <- crop(x, cropfrom)
    t00 <- ifel(t00>0,NA,NA)
    t0 <- merge(t00,t0, na.rm=FALSE); names(t0) <- 't0'
  }
  if(!is.null(elev)){
    e0 <- project(elev, t0) |> crop(cropto)
    u <- terra::linearUnits(t0)
    u <- ifelse(u == 0, 10000000/90, u)
    rs <- (res(e0)*u)[1]
    emd <- focalmed(e0, rs*10)

    erel <- e0-emd
    erel <- ifel(erel > 0,erel,0)^0.5 - ifel(-erel > 0,-erel,0)^0.5; names(erel)<-'erel'
    rss <- c(t0, e0, erel, xy0)
  }else{
    rss <- c(t0, xy0)

  }





rssnames <- names(rss)

depvar <- rssnames[1]
covars1 <- rssnames[2:5]
covars2 <- rssnames[2:length(rssnames)]
wts<- c(1,1,c(1:((rotations+1)*2))*0+1/(rotations+1))
covars3 <- c(rssnames[2:5],'resids')
df0 <- terra::spatSample(rss, 15000) |> subset(!is.na(t0))
dfsp <- st_as_sf(df0, coords = c(x='lon', y='lat'), crs=st_crs(t1))

# f.glm <- as.formula(paste(paste(depvar,paste(paste(covars1, "+ I(",covars1,"^2)",collapse = " + ", sep = ""),""), sep = " ~ ")
# ))
f.glm <- as.formula(paste(paste(depvar,paste(paste(covars1, collapse = " + ", sep = ""),""), sep = " ~ ")
))

f.rf <- as.formula(paste(paste("resids",paste(paste(covars2, collapse = " + ", sep = ""),""), sep = " ~ ")
))
# f.rf <- as.formula(paste(paste(depvar,paste(paste(covars2, collapse = " + ", sep = ""),""), sep = " ~ ")
# ))

# f.glm2 <- as.formula(paste(paste(depvar,paste(paste(covars3, "+ I(",covars3,"^2)",collapse = " + ", sep = ""),""), sep = " ~ ")
# ))
f.glm2 <- as.formula(paste(paste(depvar,paste(paste(covars3, collapse = " + ", sep = ""),""), sep = " ~ ")
))


gm <- glm(f.glm,
          family='gaussian',
          data=df0)

df0 <- df0 |> mutate(pred = predict(gm, df0), resids = t0-pred)
summary(gm)
rf <- ranger(f.rf,
             split.select.weights=wts,
             data=df0)
resids <- terra::predict(rss, rf)
resids <- focalmed(resids, rs*3); names(resids) <- 'resids'


rss1 <- c(rss, resids)

df1 <- terra::spatSample(rss1, 15000) |> subset(!is.na(t0))


gm2 <- glm(f.glm2,
          family='gaussian',
          data=df1)

pred <- terra::predict(rss1, gm2) ; names(pred) <- 'pred'}
return(pred)

plot(pred)
