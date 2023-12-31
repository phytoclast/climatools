#' @title Estimate extreme annual low temperature
#'
#' @param Tcl Mean daily low temperature of coldest month (degrees Celsius)
#' @param Lat Latitude (-90 - 90)
#' @param Lon Longitude (-180 - 180)
#' @param Elev Elevation above sea level (m).
#'
#' @return Extreme annual low temperature (degrees Celsius) as used for plant hardiness zones.
#' Estimated with a linear regression against previously published hardiness zone maps (Magarey, Borchert, & Schlegel, 2008; Daly et al., 2012) and gridded mean daily low temperature layers (Hijmans and Cameron, 2005).
#' @references Daly, C., Widrlechner, M.P., Halbleib, M.D., Smith, J.I. and Gibson, W.P., 2012. Development of a new USDA plant hardiness zone map for the United States. Journal of Applied Meteorology and Climatology, 51(2), pp.242-264.
#' @references Hijmans, R.J., S.E. Cameron, J.L. Parra, P.G. Jones and A. Jarvis, 2005. Very high resolution interpolated climate surfaces for global land areas. International Journal of Climatology 25: 1965-1978. (http://www.worldclim.org/version1)
#' @references Magarey, R.D., Borchert, D.M. and Schlegel, J.W., 2008. Global plant hardiness zones for phytosanitary risk analysis. Scientia Agricola, 65(SPE), pp.54-59.
#'
#' @export
#'
#' @examples XtremLow(15, 24, -82, 0)
XtremLow <- function(Tcl, Lat, Lon, Elev){
  Tcl <- min(Tcl)
  pacificsouth <- 1/((((Lat - -22.7)/13)^2 + ((Lon - -82.3)/14)^2)^2+1)
  amazon2 <- 1/((((Lat - -10.2)/5)^2 + ((Lon - -59.9)/10)^2)^2+1)
  amazon1 <- 1/((((Lat - -2.8)/14)^2 + ((Lon - -61.3)/19)^2)^2+1)
  pacificcent <- 1/((((Lat - 4.1)/21)^2 + ((Lon - -122.4)/41)^2)^2+1)
  mexico <- 1/((((Lat - 26)/6)^2 + ((Lon - -98.4)/12)^2)^2+1)
  florida <- 1/((((Lat - 27.5)/4)^2 + ((Lon - -81.1)/8)^2)^2+1)
  pacificnorth <- 1/((((Lat - 32.9)/26)^2 + ((Lon - -145)/27)^2)^2+1)
  oklahoma <- 1/((((Lat - 33.6)/4)^2 + ((Lon - -98.4)/8)^2)^2+1)
  arizona <- 1/((((Lat - 34)/12)^2 + ((Lon - -113.1)/8)^2)^2+1)
  atlantic <- 1/((((Lat - 34)/15)^2 + ((Lon - -60.7)/19)^2)^2+1)
  himalayas <- 1/((((Lat - 35.3)/6)^2 + ((Lon - 91.3)/13)^2)^2+1)
  kentucky <- 1/((((Lat - 38.5)/3)^2 + ((Lon - -87.6)/9)^2)^2+1)
  detroit <- 1/((((Lat - 41.8)/3)^2 + ((Lon - -82.6)/4)^2)^2+1)
  ontario <- 1/((((Lat - 44.6)/2)^2 + ((Lon - -79.2)/6)^2)^2+1)
  montana <- 1/((((Lat - 45.4)/5)^2 + ((Lon - -111.8)/10)^2)^2+1)
  minn <- 1/((((Lat - 47.6)/6)^2 + ((Lon - -92.6)/12)^2)^2+1)
  hudson <- 1/((((Lat - 60)/7)^2 + ((Lon - -87)/34)^2)^2+1)
  siberia <- 1/((((Lat - 61.2)/20)^2 + ((Lon - 105.7)/39)^2)^2+1)
  california <- 1/((((Lat - 34.8)/9)^2 + ((Lon - -128.2)/9)^2)^2+1)
  washington <- 1/((((Lat - 46)/5)^2 + ((Lon - -126.6)/5)^2)^2+1)
  colorado <- 1/((((Lat - 38.3)/2)^2 + ((Lon - -108.8)/3)^2)^2+1)
  hawaii <- 1/((((Lat - 21.3)/7)^2 + ((Lon - -157.5)/11)^2)^2+1)
  chess <- 1/((((Lat - 37)/3)^2 + ((Lon - -74)/3)^2)^2+1)

  Tclx<-	-9.171	+
    Tcl *	1.202	+
    Lat *	-0.04149	+
    Elev *	0.0008691	+
    Lat * Elev *	-0.00002455	+
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
  return(Tclx)}


#' @title Growing season temperature
#'
#' @param t Vector of monthly mean temperatures for 12 months (degrees Celsius)
#'
#' @return Mean positive temperature of warmest 6 consecutive months (degrees Celsius; all monthly temperatures below zero counted as zero)
#' @export
#'
#' @examples t <- generateTemp(-5,21)
#' GrowTemp(t)
GrowTemp <- function(t){
  bt <- ifelse(t > 0,t,0)
  tg <- pmax(mean(bt[c(1:4,11:12)]), mean(bt[c(5:10)]))
  return(tg)}

#' @title Warmest Trimester
#'
#' @param t Vector of monthly mean temperatures for 12 months (degrees Celsius)
#'
#' @return Mean positive temperature of warmest 3 consecutive months (degrees Celsius; all monthly temperatures below zero counted as zero)
#' @export
#'
#' @examples t <- generateTemp(-5,21)
#' GrowTemp(t)
SummerTemp <- function(t){
  bt <- ifelse(t > 0,t,0)
  bt <- c(bt,bt[1:2])
  for(i in 1:12){
    bt0<- mean(bt[i:(i+2)])
    if(i==1){bt3max<-bt0}else{bt3max<-pmax(bt3max,bt0)}
  }
  return(bt3max)}


#' @title Mean annual biotemperature
#'
#' @param t Vector of monthly mean temperatures for 12 months (degrees Celsius)
#'
#' @return Mean annual temperature between 0 and 30 (degrees Celsius; all monthly temperatures below zero counted as zero, and those above 30 limited to 30, prior to averaging)
#' @references Holdridge, L.R. and Tosi Jr, J.A., 1967. Tropical Science Center: San jose, Costa Rica. Life Zone Ecology; Tropical Science Center: San Jose, Costa Rica.
#' @export
#'
#' @examples t <- generateTemp(-5,21)
#' Biotemperature(t)
Biotemperature <- function(t){
  bt <- mean(ifelse(t > 0,ifelse(t>30,30,t),0))
  return(bt)}

#' @title Get peak 3-month evapotranspiration
#' @description
#' Function calculates p3AET, the peak consecutive 3-month sum of the lesser among monthly precipitation and evapotranspiration. Ignoring soil water storage, this is equivalent to the expected actual evapotranspiration of the most productive 3 months. This is a measure of the degree to which precipitation is associated with warm temperatures. A Mediterranean climate with the same temperatures and precipitation as a subtropical savanna will tend to have a lower value than the savanna. Higher value would tend to associate with a greater proportion of precipitation as thunderstorms.
#'
#'
#' @param p vector of 12 monthly precipitation (mm)
#' @param e vector of 12 monthly potential evapotranspiration (mm)
#'
#' @return p3AET (mm)
#' @export
#'
#' @examples mon <- 1:12
#' t <- generateTemp(-5,21)
#' #Wet summer
#' p <- generatePpt(700/12, 160, 5, 6, 12)
#' e <- GetPET(mon,43,t+5,t-5,p)
#' Getp3AET(p,e)
#' #Wet winter
#' p <- generatePpt(700/12, 160, 5, 12, 6)
#' e <- GetPET(mon,43,t+5,t-5,p)
#' Getp3AET(p,e)
Getp3AET <- function(p,e){
  paet <- pmin(p,e)
  paet1 <- c(paet[11:12],paet,paet[1:2])
  for(i in 1:12){
    paet0 <- sum(paet1[(i):(i+2)])
  if(i==1){p3AET=paet0}else{
    p3AET=pmax(p3AET, paet0)
  }}
  return(p3AET)}


#' Highest/Lowest n month period temperature/precipitation
#'
#' @param t vector of 12 monthly mean temperatures
#' @param p vector of 12 monthly mean precipitation
#' @param n number of months in period
#' @param max_var Use maximum if TRUE (minimum if FALSE)
#' @param v_input Assess which parameter to rate extreme for period (t = temperature; p = precipitation)
#' @param v_return Return which parameter for the extreme period if TRUE (t = temperature; p = precipitation; b = positive biotemperature)
#' @param bt Temperature threshold for biotemperature
#'
#' @return A single value for the temperature/precipitation of the warmest/coldest/wettest/driest n-month period.
#' @export
#'
#' @examples t=generateTemp(5,22)
#' p=generatePpt(100,50,150)
#' #precipitation of warmest 6 months
#' warm6ppt(t, p, n=6, max_var=TRUE, v_input='t', v_return='p')
#' #temperature of warmest 6 months
#' warm6ppt(t, p, n=6, max_var=TRUE, v_input='t', v_return='t')
#' #temperature of driest 3 months
#' warm6ppt(t, p, n=3, max_var=FALSE, v_input='p', v_return='t')
#' #precipitation of driest 3 months
#' warm6ppt(t, p, n=3, max_var=FALSE, v_input='p', v_return='p')
#' #temperature of wettest 3 months
#' warm6ppt(t, p, n=3, max_var=TRUE, v_input='p', v_return='t')
#' #positive biotemperature of warmest 6 months (default threshold of 0 degrees Celsius)
#' warm6ppt(t, p, n=6, max_var=TRUE, v_input='t', v_return='b')
#' #positive biotemperature of warmest 6 months, with a threshold of 5 degrees Celsius
#' warm6ppt(t, p, n=6, max_var=TRUE, v_input='t', v_return='b', bt = 5)
#' #positive biotemperature of 12 months
#' warm6ppt(t, p, n=12, max_var=TRUE, v_input='t', v_return='b')
#' #positive biotemperature of 12 months, with a threshold of 5 degrees Celsius
#' warm6ppt(t, p, n=12, max_var=TRUE, v_input='t', v_return='b', bt = 5)
#'
warm6ppt <- function(t,p,n=6,max_var=TRUE,v_input=c('t','p'),v_return=c('t','p','b'),bt=0){
  #v_input=c('t','p')[1];v_return=c('t','p','b')[3]
  ctab <- data.frame(t=t,p=p)
  if(v_input %in% 't'){ctab$v = ctab$t}else{ctab$v = ctab$p}
if(v_return %in% 't'){ctab$r = ctab$t}else if(v_return %in% 'b'){
  ctab$r = ifelse(ctab$t > bt, ctab$t - bt,0)
}else{ctab$r = ctab$p}
  ctab$v6 <- NA_real_
  ctab$r6 <- NA_real_
  ctab <- rbind(ctab,ctab)
  for(i in 1:12){#i=1
    ctab[i,]$v6 <- sum(ctab[i:(i+n-1),]$v)
    ctab[i,]$r6 <- sum(ctab[i:(i+n-1),]$r)
  }
  ctab <- ctab[1:12,]
  if(max_var){extv <- max(ctab$v6)}else{extv <- min(ctab$v6)}
  extr <- subset(ctab, v6 %in% extv)$r6[1]
  if(v_return %in% c('t','b')){extr <- extr/n}
  return(extr)
}
