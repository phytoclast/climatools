library(terra)
library(climatools)
altlayer <- rast(system.file("extdata", "altlayer.tif", package="climatools") )
pts.p <- read.csv(system.file("extdata", "pts.p.csv", package="climatools"))
pts.t <- read.csv(system.file("extdata", "pts.t.csv", package="climatools"))
pts <- pts.p[,c('lon','lat','p01')]
pts[,3] <- log10(pts[,3] + 1)
p01 <- toclimrast(pts, altlayer, covrange = 500, segx = 10, segy = 10)
pts <- pts.p[,c('lon','lat','p07')]
pts[,3] <- log10(pts[,3] + 1)
p07 <- toclimrast(pts, altlayer, covrange = 500, segx = 10, segy = 10)
pts <- pts.t[,c('lon','lat','t01','elev')]
t01 <- toclimrast(pts, altlayer, covrange = 500, segx = 10, segy = 10)
pts <- pts.t[,c('lon','lat','t07','elev')]
t07 <- toclimrast(pts, altlayer, covrange = 500, segx = 10, segy = 10)
climrast <- c(p01,p07,t01,t07);names(climrast) <- c('p01','p07','t01','t07')
plot(climrast)
grd <- c(altlayer,climrast)
pts <- pts.t |> mutate(pos = ifelse((lon > -110 & lon < -100 & lat > 30 & lat < 45 &
                                      elev > 1000 & elev < 1500) | (lat < 40 & elev <100), 1,0))

vts <- vect(pts[pts$pos==1,],geom=c("lon", "lat"),crs=crs('epsg:4326'))
plot(t01);points(vts)

pts <- pts[,c('lon','lat','pos')]
altlayer=grd; cropto=NULL; covrange=0; minrow=50; segx=10; segy=10;
cropbuffer=5; randforest = TRUE; smoothresiduals = TRUE; refit=FALSE

geoglm <- function(pts, altlayer, cropto=NULL, covrange=0, minrow=50, segx=5, segy=5,
                       cropbuffer=5, randforest = TRUE, smoothresiduals = TRUE, refit=FALSE){
  #crop full extent if null
  if(is.null(cropto)){cropto <- terra::ext(altlayer)
  cropbuffer=0}
  #ensure that input has only 3 columns
  pts <- as.data.frame(pts)
  if(ncol(pts)>3){
    pts <- pts[,1:4]
    names(pts) <- c('x','y','z','original')
  }else{
    pts <- pts[,1:3]
    names(pts) <- c('x','y','z')
  }
  #crop raster to new extent with buffer
  cropto0 <- cropto + c(-cropbuffer,cropbuffer,-cropbuffer,cropbuffer)
  #crop point data extent and convert to terra spatial vector.
  pts <- pts |> subset(x >= cropto0[1] & x <= cropto0[2] &
                         y >= cropto0[3] & y <= cropto0[4])
  vts <- vect(pts, geom=c("x", "y"),crs=crs('epsg:4326'))

  grd <- crop(altlayer, cropto0)
  #add raster with rotated xy coordinates
  xy0 <- climatools::makexyrast(grd[[1]],6)
  grdall <- c(grd, xy0)
  #create lower resolution raster to model covariates and residuals
  grdall.1 <- aggregate(grdall,fact=3,fun="mean")
  #get raster units to ensure that focal analyses neighborhoods are consistent
  u <- terra::linearUnits(grd)
  u <- ifelse(u == 0, 10000000/90, u)
  rs <- (terra::res(grd)*u)[1]

  #extract rasters to points
  vts <- project(vts,altlayer)
  vtsgrd <- terra::extract(grdall,vts)
  pts <- cbind(pts,vtsgrd)
  #if original station elevation provided, transfer to first covariate
  if('original' %in% colnames(pts)){
    firstcov <- names(grd)[1]
    pts[,firstcov] <- pts[,'original']}

  #build formulas
  depvar <- names(pts)[3]
  covars1 <- c(names(grd),names(xy0)[1:2])
  covars2 <- c(names(grd),names(xy0))
  # f.glm <- stats::as.formula(paste(paste(depvar,paste(paste("poly(",covars1,",2)", collapse = " + ", sep = ""),""), sep = " ~ ")
  # ))

  #prepare regression loops
  exfactors <- c(0,0.5,1,2,5,10)
  spanx <- (cropto[2]-cropto[1])/segx
  spany <- (cropto[4]-cropto[3])/segy
  pts$inner <- NA
  pts$outer <- NA
  pts$coeffs0 <- NA
  pts$erange <- NA
  #loop through x and y segments
  for(i.x in 1:segx){
    for(i.y in 1:segy){
      #set status for this segment until acceptable regression model
      success <- FALSE
      #gradually expand size of analysis area
      for(i.f in 1:length(exfactors)){
        #i.x=3;i.y=4;i.f=1
        if(!success){
          exfact <- exfactors[i.f]
          addtoy <- spany*exfact
          addtox <- spanx*exfact
          #establish inner and outer points; inner points of segment carry the values of covariates; outer values are involved in the building the models
          crop0 <- c(cropto[1]+(i.x-1)*spanx,cropto[1]+i.x*spanx,
                     cropto[3]+(i.y-1)*spany,cropto[3]+i.y*spany)
          pts <- pts |> mutate(inner = ifelse(x >= crop0[1] & x <= crop0[2] &
                                                y >= crop0[3] & y <= crop0[4], 1, 0),
                               outer = ifelse(x >= (crop0[1]-addtox) & x <= (crop0[2]+addtox) &
                                                y >= (crop0[3]-addtoy) & y <= (crop0[4]+addtoy), 1, 0))
          pts.i <- pts |> subset(outer ==1)
          #move on if segment is has too few rows to build model
          if(nrow(pts.i) > minrow){
            #move on (expand extent) if segment points do not have sufficient range in first variable (e.g. elevation) to build accurate model
            erange0 <- nrow(pts[pts.i$z %in% 1,])
            if(erange0 >= 10){
              #model segment of points and feed coefficients into points dataset

              cov.unique <- apply(pts.i[,covars1], MARGIN=2, FUN=function(x){length(unique(x))})
              usecovs1 <- names(cov.unique[cov.unique >= 2])
              usecovs2 <- names(cov.unique[cov.unique > 3])

              f.glm <- stats::as.formula(paste(paste(depvar,paste(paste(usecovs1,"+","I(",usecovs2,"^2)", collapse = " + ", sep = ""),""), sep = " ~ ")
              ))
              gm <- stats::glm(f.glm,
                               family='binomial',
                               data=pts.i)
              summary(gm)
              cofs <- list(gm$coefficients)
              pts <- pts |> mutate(coeffs0 = ifelse(inner %in% 1, cofs,coeffs0))
              success <- TRUE}
          }}}
    }}



  #extract coefficients from points and convert to rasters using either randomforest model or interpolation

 allcn <- c("(Intercept)", covars1, paste0("I(",covars1,"^2)"))
 allcndf <- data.frame(coffs=0,rownams=allcn)
 cflist <- matrix(nrow = nrow(pts), ncol = length(covars1)*2+1)
 cflist <- as.data.frame(cflist)
 colnames(cflist) <- allcn
 for(i.cv in 1:nrow(pts)){
  l <- pts$coeffs0[i.cv]
  df <- as.data.frame(l)
  colnames(df) <- 'coffs'
  df$rownams <- rownames(df)
  df <- rbind(df,allcndf)
  df <- aggregate(coffs ~ rownams, data=df, FUN=sum)
  dfrn <- df$rownams
  df <- t(df[,2])
  colnames(df) <- dfrn
  df <- df[,allcn,drop = FALSE]
  cflist[i.cv,] <- df}
  nc <- ncol(cflist)
  grdall.0 <- rast(resolution=res(grdall.1), crs=crs(grdall.1), extent=ext(grdall.1), nlyrs=nc)
  for(i in 1:nc){
    pts$coeffs <- cflist[,i]
    if(randforest){
      f.rf <- stats::as.formula(paste(paste("coeffs",paste(paste(covars2, collapse = " + ", sep = ""),""), sep = " ~ ")
      ))
      rf <- ranger::ranger(f.rf,
                           # split.select.weights=wts,
                           num.trees = 100,
                           data=pts[!is.na(pts$coeffs),])

      cofffs <- terra::predict(grdall.1, rf)
    }else{
      xyz <- pts[,c('x','y','coeffs')] |> subset(!is.na(coeffs))
      gs <- gstat::gstat(formula=coeffs~1, locations=~x+y, data=xyz, nmax=32, set=list(idp = 2))
      cofffs <- interpolate(grdall.1, gs, debug.level=0)[[1]]
    }
    cofffs <- focalmed(cofffs, segy*u/3);
    names(cofffs)  <- paste0("coef.",i)
    grdall.0[[i]] <- cofffs
  }
  grdall.0 <- project(grdall.0, grdall)

  #extract coefficients to points
  pts2 <- pts |> cbind(extract(grdall.0,vts))
  grdall2 <- c(grdall, grdall.0)
  #create formula with covariates and coefficients
  # ii= 1:nc
  # odds <- (ii)/2!=floor((ii)/2)
  # cofodds <- ii[odds][-1]
  # cofevens <- ii[!odds]
  covarc2 <- names(grdall.0)[((nc-1)/2+2):nc]
  covarc1 <- names(grdall.0)[2:((nc-1)/2+1)]
  intcp <- names(grdall.0)[1]





  #formula changed to only have interaction and not coefficient as its own covariate
  f.glm2 <- stats::as.formula(paste(depvar,paste(intcp, paste("poly(",covars1,",2)",
                                                              "+",covars1,":",covarc1,
                                                              "+","I(",covars1,"^2):",covarc2, collapse = " + ", sep = ""), sep = " + "), sep = " ~ "))
  f.rf <- stats::as.formula(paste(depvar,paste(intcp, paste("poly(",covars1,",2)", collapse = " + ", sep = ""), sep = " + "), sep = " ~ "))

  #linear model with new formula
  gm2 <- stats::glm(f.glm2,
                    family='binomial',
                    data=pts2)
  # summary(gm2)
   1-gm2$deviance/gm2$null.deviance

   #use model to generate prediction layer
  pred <- terra::predict(grdall2, gm2, type='response')
plot(pred, col=map.pal('bcyr'));points(vts[vts$z==1])

f.rf <- stats::as.formula(paste(paste("z",paste(paste(covars2, collapse = " + ", sep = ""),""), sep = " ~ ")
))
rf <- ranger::ranger(f.rf,
                     # split.select.weights=wts,
                     #num.trees = 1500,
                     data=pts2)
pred <- terra::predict(grdall2, rf)
plot(pred, col=map.pal('bcyr'));points(vts[vts$z==1])
1-rf$prediction.error
  return(pred)
}
