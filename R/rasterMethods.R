#Raster methods for climatic functions

GetSolarRad.rast <- function(DayNumber, Lat){
  declination <- GetDcl(DayNumber)

  hs <- acos(min(max(-tan(Lat/360*2*3.141592) * tan(declination),-1),1))
  Ra <- 117.5 * (hs*sin(Lat/360*2*3.141592)*sin(declination) +
                   cos(Lat/360*2*3.141592)*cos(declination)*sin(hs)) / 3.141592
  return(Ra)
}
GetPET.rast <- function(mon, lat, th, tl, p){
  Ra=GetSolarRad.rast(GetDayNumber(mon),lat)
  e =  0.85*GetTransGrow.rast(th, tl)*GetPETdaily.rast(Ra, th, tl, p)*GetNumberDay(mon) #0.85829 is crop coefficient to make comparable to Thornthwaite. Multiplied by growing season to zero out transpiration portion of evapotranspiration in freezing weather. Per day rate need to be multiplied by days per month to get monthly rate.
  mons <- c('e01','e02','e03','e04','e05','e06','e07','e08','e09','e10','e11','e12')
  names(e) <- mons[mon]
  return(e)}

GetTransGrow.rast <- function(th, tl) {#Adjust to reduction in transpiration due to cold, with evaporation only outside growing season
  ts = 0.8 #assumed T/ET ratio during growing season
  tw = 0 #assumed T/ET ratio during freezing season
  t <- (th+tl)/2
  tr <- 10 #generally as mean temperatures get below 10 transpiration shuts down, regardless of warm daytime temperatures
  G0 <- (t-0)/(tr)
  G1 <- min(max(G0,0),1) #generally as mean temperatures get below 5 transpiration shuts down, regardless of warm daytime temperatures
  evmin = (tw)+(1-ts)
  G = G1*(1-evmin)+evmin
  return(G)}

GetPETdaily.rast <- function(Ra, th, tl, p){
  Vpmax = 0.6108*exp(17.27*th/(th+237.3)) #saturation vapor pressure kPa
  Vpmin = 0.6108*exp(17.27*tl/(tl+237.3)) #saturation vapor pressure kPa
  logp <- log(p+1)
  e0 <- Ra*0.0508780  +
    Vpmax*0.7893714  +
    Vpmin*-0.5589255  +
    logp*-0.1309403  +
    Ra*Vpmax*0.0049383
  e <- max(e0,0)
  return(e)}

GetPET.block <- function(mon, block, lat='lat', th.jan='th01', tl.jan='tl01', p.jan='p01'){
  th.ind = which(names(block) %in% th.jan)
  tl.ind = which(names(block) %in% tl.jan)
  p.ind = which(names(block) %in% p.jan)
  lat.ind = which(names(block) %in% lat)
  th <- block[,,(th.ind+mon-1),drop=FALSE]
  tl <- block[,,(tl.ind+mon-1),drop=FALSE]
  p  <- block[,,(p.ind+mon-1),drop=FALSE]
  Lat  <- block[,,(lat.ind),drop=FALSE]
  e <- GetPET.rast(mon=mon, lat = Lat, th= th, tl= tl, p= p)
  return(e)
}



meanT.rast <- function(block, th.jan='th01', tl.jan='tl01'){
  th.ind = which(names(block) %in% th.jan)
  tl.ind = which(names(block) %in% tl.jan)
  th <- block[,,th.ind:(th.ind+11),drop=FALSE]
  tl <- block[,,tl.ind:(tl.ind+11),drop=FALSE]

  t <- (th+tl)/2

  names(t) <- c('t01','t02','t03','t04','t05','t06','t07','t08','t09','t10','t11','t12')
  return(t)}


GrowTemp.rast <- function(block, t.jan='t01'){
  t.ind = which(names(block) %in% t.jan)
  t.val <- block[,,t.ind:(t.ind+11),drop=FALSE]
  t.val <- ifel(t.val>0,t.val,0)
  t.val <- max(mean(t.val[,,c(11,12,1,2,3,4),drop=FALSE]),mean(t.val[,,c(5,6,7,8,9,10),drop=FALSE]))
  names(t.val) <- 'Tg'
  return(t.val)
  }

#summarizes a climate statistic from a raster stack of monthly stats consecutively arranged by month; user identifies name of first month of that statistic and name of the function to summarize it.
ApplyClim.rast <- function(block, jan='p01',mons=c(1,2,3,4,5,6,7,8,9,10,11,12),  fun='sum', name = NULL){
  ind = which(names(block) %in% jan)
  x <- block[,,ind:(ind+11),drop=FALSE]
  x <- x[,,mons, drop=FALSE]
  x <- terra::app(x, fun=fun)
  if(!is.null(name)){
    names(x) <- name
  }
  return(x)
}

#Calculates the maximum possible evapotranspiration based only on precipitation for each month, without monthly carry over of soil water storage.
AET.rast <- function(block, jan.p='p01',jan.e='e01'){
  p.ind = which(names(block) %in% jan.p)
  e.ind = which(names(block) %in% jan.e)
  p <- block[,,p.ind:(p.ind+11),drop=FALSE]
  e <- block[,,e.ind:(e.ind+11),drop=FALSE]
  a <- e
  for(i in 1:12){
  a[[i]] <- min(e[[i]],p[[i]])}
  x <- terra::app(a, fun='sum')
  names(x) <- 'a'
  return(x)
}

AET.rast.monthly <- function(block, jan.p='p01',jan.e='e01'){
  month <- c('01','02','03','04','05','06','07','08','09','10','11','12')
  p.ind = which(names(block) %in% jan.p)
  e.ind = which(names(block) %in% jan.e)
  p <- block[,,p.ind:(p.ind+11),drop=FALSE]
  e <- block[,,e.ind:(e.ind+11),drop=FALSE]
  a <- e
  for(i in 1:12){
    a[[i]] <- min(e[[i]],p[[i]])
    names(a[[i]]) <- paste0('a',month[i])
  }
  return(a)
}

AET.rast.max <- function(block, jan.a='a01', nmonth = 3){
  a.ind = which(names(block) %in% jan.a)
  a <- block[,,a.ind:(a.ind+11),drop=FALSE]
  indices <- c(1:12,1:12)
  amx <- a
  for(i in 1:12){
    thismonth <- indices[i:(i+nmonth-1)]
    amx[[i]] <- sum(a[[thismonth]])
  }
  x <- max(amx)
  names(x) <- paste0('max',nmonth,'aet')
  return(x)
}



XtremLow.rast <- function(Tcl, elev){
  Tcl <- min(Tcl)
  xyr <- rast(nrows=180, ncols=360, nlyrs=1, crs=crs('EPSG:4326'), extent=ext(-180,180,-90,90), vals=1)
  xy <- terra::as.data.frame(x=xyr, xy=TRUE)
  lon <- rast(x=xy[,c('x','y','x')], type="xyz", crs=crs(xyr), extent=ext(xyr))
  lat <- rast(x=xy[,c('x','y','y')], type="xyz", crs=crs(xyr), extent=ext(xyr))
  names(lon) <- 'lon'
  names(lat) <- 'lat'

  pacificsouth <- 1/((((lat - -22.7)/13)^2 + ((lon - -82.3)/14)^2)^2+1)
  amazon2 <- 1/((((lat - -10.2)/5)^2 + ((lon - -59.9)/10)^2)^2+1)
  amazon1 <- 1/((((lat - -2.8)/14)^2 + ((lon - -61.3)/19)^2)^2+1)
  pacificcent <- 1/((((lat - 4.1)/21)^2 + ((lon - -122.4)/41)^2)^2+1)
  mexico <- 1/((((lat - 26)/6)^2 + ((lon - -98.4)/12)^2)^2+1)
  florida <- 1/((((lat - 27.5)/4)^2 + ((lon - -81.1)/8)^2)^2+1)
  pacificnorth <- 1/((((lat - 32.9)/26)^2 + ((lon - -145)/27)^2)^2+1)
  oklahoma <- 1/((((lat - 33.6)/4)^2 + ((lon - -98.4)/8)^2)^2+1)
  arizona <- 1/((((lat - 34)/12)^2 + ((lon - -113.1)/8)^2)^2+1)
  atlantic <- 1/((((lat - 34)/15)^2 + ((lon - -60.7)/19)^2)^2+1)
  himalayas <- 1/((((lat - 35.3)/6)^2 + ((lon - 91.3)/13)^2)^2+1)
  kentucky <- 1/((((lat - 38.5)/3)^2 + ((lon - -87.6)/9)^2)^2+1)
  detroit <- 1/((((lat - 41.8)/3)^2 + ((lon - -82.6)/4)^2)^2+1)
  ontario <- 1/((((lat - 44.6)/2)^2 + ((lon - -79.2)/6)^2)^2+1)
  montana <- 1/((((lat - 45.4)/5)^2 + ((lon - -111.8)/10)^2)^2+1)
  minn <- 1/((((lat - 47.6)/6)^2 + ((lon - -92.6)/12)^2)^2+1)
  hudson <- 1/((((lat - 60)/7)^2 + ((lon - -87)/34)^2)^2+1)
  siberia <- 1/((((lat - 61.2)/20)^2 + ((lon - 105.7)/39)^2)^2+1)
  california <- 1/((((lat - 34.8)/9)^2 + ((lon - -128.2)/9)^2)^2+1)
  washington <- 1/((((lat - 46)/5)^2 + ((lon - -126.6)/5)^2)^2+1)
  colorado <- 1/((((lat - 38.3)/2)^2 + ((lon - -108.8)/3)^2)^2+1)
  hawaii <- 1/((((lat - 21.3)/7)^2 + ((lon - -157.5)/11)^2)^2+1)
  chess <- 1/((((lat - 37)/3)^2 + ((lon - -74)/3)^2)^2+1)

  latlon <-	-9.171	+
    lat *	-0.04149	+
    pacificsouth *	-1.792	+
    amazon2 *	2.573	+
    amazon1 *	-1.014	+
    pacificcent *	-0.749	+
    mexico *	-0.8227	+
    florida *	-3.557	+
    pacificnorth *	-1.246	+
    oklahoma *	0.1758	+
    arizona *	2.605	+
    chess *	0.8347	+
    atlantic *	0.2967	+
    himalayas *	-1.814	+
    kentucky *	-2.644	+
    detroit *	0	+
    ontario *	-2.314	+
    montana *	-4.415	+
    minn *	1.136	+
    hudson *	-5.154	+
    siberia *	-3.797	+
    california *	4.48	+
    washington *	3.597	+
    colorado *	1.458	+
    hawaii *	6.673

  lat <- project(lat, elev)

  latlon <- project(Tclx0, elev)

  Tclx <-	Tcl *	1.202	+
    latlon +
    lat *	-0.04149	+
    elev *	0.0008691	+
    lat * elev *	-0.00002455
  names(Tclx) <- 'Tclx'
  return(Tclx)}

