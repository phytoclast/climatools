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
  t <- block[,,t.ind:(t.ind+11),drop=FALSE]
  t <- ifel(t>0,t,0)
  t <- max(mean(t[,,c(11,12,1,2,3,4),drop=FALSE]),mean(t[,,c(5,6,7,8,9,10),drop=FALSE]))
  names(t) <- 'Tg'
  return(t)
}

Warm.rast <- function(block, t.jan='t01'){
  t.ind = which(names(block) %in% t.jan)
  t <- block[,,t.ind:(t.ind+11),drop=FALSE]
  t <- max(t)
  names(t) <- 'Tw'
  return(t)
}
WarmHigh.rast <- function(block, th.jan='th01'){
  t.ind = which(names(block) %in% th.jan)
  t <- block[,,t.ind:(t.ind+11),drop=FALSE]
  t <- max(t)
  names(t) <- 'Twh'
  return(t)
}
Cold.rast <- function(block, t.jan='t01'){
  t.ind = which(names(block) %in% t.jan)
  t <- block[,,t.ind:(t.ind+11),drop=FALSE]
  t <- min(t)
  names(t) <- 'Tc'
  return(t)
}
ColdLow.rast <- function(block, tl.jan='tl01'){
  t.ind = which(names(block) %in% tl.jan)
  t <- block[,,t.ind:(t.ind+11),drop=FALSE]
  t <- min(t)
  names(t) <- 'Tcl'
  return(t)
}