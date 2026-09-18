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
x <- t1#thisgrid#; t1<-logp(t1)
elev <- dem
# cropto <- c(-90, -75, 40, 55)
# cropfrom <- c(-82,-80, 45.5,49)
cropto <- c(-152,-149.5, 62.6, 63.2)
cropfrom <- c(-151.5,-150, 62.8, 63)


plot(x)
x2 <- refit(x, elev=dem, cropto=cropto, rotations=6)
plot(x2)






refit <- function(x, elev=NULL, cropto=NULL, cropfrom=NULL, rotations=0){
  #If cropping extent provided
  if(!is.null(cropto)){
    t0 <- crop(x, cropto);
  }
  #Standardize name for extracting to a data frame to be used by formula
  names(t0) <- 't0'

  #dummy values for having no data for these layers
  erel <- NULL
  wt1 <- NULL; wt2 <- NULL

  #Create alternative rotated XY coordinates to make a smoother random forest model.
  xy0 <- makexyrast(t0,rotations)

  #Omit problematic data to patch with model using less problematic data from adjacent area.
  if(!is.null(cropfrom)){
    t00 <- crop(x, cropfrom)
    t00 <- ifel(t00>0,NA,NA)
    t0 <- merge(t00,t0, na.rm=FALSE); names(t0) <- 't0'
  }

  #If elevation data is used, create auxiliary layer for relative elevation which captures local inversions.
  if(!is.null(elev)){
    e0 <- project(elev, t0) |> crop(cropto)
    u <- terra::linearUnits(t0)
    u <- ifelse(u == 0, 10000000/90, u)
    rs <- (res(e0)*u)[1]
    emd <- focalmed(e0, rs*10)

    erel <- e0-emd
    erel <- ifel(erel > 0,erel,0)^0.5 - ifel(-erel > 0,-erel,0)^0.5; names(erel)<-'erel'
    rss <- c(t0, e0, erel, xy0)
    wt1 <- 1; wt2 <- 1
  }else{
    #Use only xy data instead
    rss <- c(t0, xy0)
  }

  #Define formulas for general linear and random forest models.
  depvar <- "t0"
  covars1 <- c(names(xy0)[1:2],names(elev),names(erel))
  covars2 <- c(names(xy0),names(elev),names(erel))
  wts<- c(c(1:((rotations+1)*2))*0+1/(rotations+1), wt1,wt2)
  covars3 <- c(covars1,'resids')

  #Sample rasters to train models (omitting missing data)
  df0 <- terra::spatSample(rss, 1500) |> subset(!is.na(t0))

  #data points converted to spatial features to test coverage
  #dfsp <- sf::st_as_sf(df0, coords = c(x='lon', y='lat'), crs=sf::st_crs(t1)); plot(vect(dfsp))

  f.glm <- as.formula(paste(paste(depvar,paste(paste(covars1, collapse = " + ", sep = ""),""), sep = " ~ ")
  ))

  f.rf <- as.formula(paste(paste("resids",paste(paste(covars2, collapse = " + ", sep = ""),""), sep = " ~ ")
  ))

  f.glm2 <- as.formula(paste(paste(depvar,paste(paste(covars3, collapse = " + ", sep = ""),""), sep = " ~ ")
  ))

  #Run initial model to get linear trends
  gm <- glm(f.glm,
            family='gaussian',
            data=df0)

  df0 <- df0 |> mutate(pred = predict(gm, df0), resids = t0-pred)
  #Create additional layer with residuals using random forest model.
  rf <- ranger(f.rf,
               split.select.weights=wts,
               #num.trees = 1500,
               data=df0)
  resids <- terra::predict(rss, rf)
  resids <- focalmed(resids, rs*3); names(resids) <- 'resids'


  rss1 <- c(rss, resids)

  df1 <- terra::spatSample(rss1, 1500) |> subset(!is.na(t0))


  gm2 <- glm(f.glm2,
             family='gaussian',
             data=df1)

  pred <- terra::predict(rss1, gm2) ; names(pred) <- 'pred'
  return(pred)}


