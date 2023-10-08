
#' Get mean day of the year number by month
#'
#' @param mon Month of the year number (1-12).
#'
#' @return Day of the year at mid month (16-350).
#' @export
#'
#' @examples GetDayNumber(6)
GetDayNumber <- function(mon){#get mean day number by month
  daynumber <- c(16.000,45.625,75.250,106.125,136.250,166.750,197.250,228.250,258.750,289.250,319.750,350.250)
  return(daynumber[mon])
}


#' Get number of days in month.
#'
#' @param mon  Month of the year number (1-12).
#'
#' @return Numbers of days in the month (for average year, 28.25-31).
#' @export
#'
#' @examples GetNumberDay(6)
GetNumberDay <- function(mon){#get number days per month
  numberday <- c(31.00,28.25,31.00,30.00,31.00,30.00,31.00,31.00,30.00,31.00,30.00,31.00)
  return(numberday[mon])
}



#' Estimate mean vapor pressure
#'
#' @param p Mean monthly precipitation (mm)
#' @param th Mean daily high temperature (degrees Celsius)
#' @param tl Mean daily low temperature (degrees Celsius)
#'
#' @return Mean vapor pressure for the day or month
#' Based on a linear regression using 10 minute WorldClim 2.0 data vapor pressure estimates.
#' @export
#'
#' @examples GetVp(100,25,15)
GetVp  <- function(p,th,tl) {#Based on linear regression using 10 minute WorldClim 2.0 data with vapor pressure estimates
  Vpmax = 0.6108*exp(17.27*th/(th+237.3)) #saturation vapor pressure kPa
  Vpmin = 0.6108*exp(17.27*tl/(tl+237.3)) #saturation vapor pressure kPa
  Vp0 <- (Vpmin*7.976e-01+
            Vpmin*log(p+1)*9.499e-02+
            Vpmin*Vpmax*-6.599e-02)
  Vp <- pmax(0,pmin(Vpmin,Vp0))
  return(Vp)}


#' Estimate solar radiation reaching the earth's surface
#'
#' @param Ra Calculated daily solar radiation hitting the earth's atmosphere
#' @param Elev Elevation above sea level (m).
#' @param th Mean daily high temperature (degrees Celsius)
#' @param tl Mean daily low temperature (degrees Celsius)
#' @param p Mean monthly precipitation (mm)
#'
#' @return Estimated solar radiation reaching the earth's surface. Based on a linear regression of 10 minute WorldClim 2.0 solar radiation estimates.
#' @export
#'
#' @examples GetSolar(GetSolarRad(GetDayNumber(6), 43), 210, 25, 15, 100)
GetSolar <- function(Ra, Elev, th, tl, p) {#Based on linear regression using 10 minute WorldClim 2.0 data with solar radiation estimates
  Rso <- (0.75+2*10^-5*Elev)*Ra
  Rs0 <- (Rso*9.521e-01+
            Rso*log(p+1)*-9.087e-02+
            Rso*tl*-3.644e-03+
            Rso*log(p+1)*th*1.335e-03)
  Rs <- pmax(0.3*Rso,pmin(Rso,Rs0))
  return(Rs)}


#' Get solar radiation reaching the top of earth's atmosphere
#'
#' @param DayNumber Day of the year (1-365)
#' @param Lat Latitude (-90 - 90)
#'
#' @return Daily solar radiation reaching the top of earth's atmosphere based on a formula used in Walter et al (2000).
#'
#' @references Walter, I.A., Allen, R.G., Elliott, R., Jensen, M.E., Itenfisu, D., Mecham, B., Howell, T.A., Snyder, R., Brown, P., Echings, S. and Spofford, T., 2000. ASCE's standardized reference evapotranspiration equation. In Watershed management and operations management 2000 (pp. 1-11).
#' @export
#'
#' @examples GetSolarRad(GetDayNumber(6), 43)
GetSolarRad <- function(DayNumber, Lat){
  declination <- GetDcl(DayNumber)

  hs <- acos(pmin(pmax(-tan(Lat/360*2*3.141592) * tan(declination),-1),1))
  Ra <- 117.5 * (hs*sin(Lat/360*2*3.141592)*sin(declination) +
                   cos(Lat/360*2*3.141592)*cos(declination)*sin(hs)) / 3.141592
  return(Ra)
}


#' Get day length
#'
#' @param DayNumber Day of the year (1-365)
#' @param Lat Latitude (-90 - 90)
#'
#' @return Number of hours of daylight in a day.
#' @export
#'
#' @examples GetDayLength(GetDayNumber(6),43)
GetDayLength<- function(DayNumber, Lat){
  declination <- GetDcl(DayNumber)

  Dl <- ifelse(Lat + declination*360/2/3.141592 > 89.16924, 24, ifelse(Lat - declination*360/2/3.141592 >= 90, 0, (atan(-((sin(-0.83/360*2*3.141592)-sin(declination)*sin(Lat/360*2*3.141592))/(cos(declination)*cos(Lat/360*2*3.141592)))/(-((sin(-0.83/360*2*3.141592)-sin(declination)*sin(Lat/360*2*3.141592))/(cos(declination)*cos(Lat/360*2*3.141592)))*((sin(-0.83/360*2*3.141592)-sin(declination)*sin(Lat/360*2*3.141592))/(cos(declination)*cos(Lat/360*2*3.141592)))+1)^0.5)+2*atan(1))/3.141592*24))
  return(Dl)}


#' Get Net Solar Radiation
#'
#' @param Ra Calculated daily solar radiation hitting the earth's atmosphere
#' @param Elev Elevation above sea level (m).
#' @param th Mean daily high temperature (degrees Celsius)
#' @param tl Mean daily low temperature (degrees Celsius)
#' @param p Mean monthly precipitation (mm)
#'
#' @return Net Daily Solar Radiation based on formulas used in Walter et al (2000).
#'
#' @references Walter, I.A., Allen, R.G., Elliott, R., Jensen, M.E., Itenfisu, D., Mecham, B., Howell, T.A., Snyder, R., Brown, P., Echings, S. and Spofford, T., 2000. ASCE's standardized reference evapotranspiration equation. In Watershed management and operations management 2000 (pp. 1-11).
#' @export
#'
#' @examples GetNetSolar(GetSolarRad(GetDayNumber(6), 43), 210, 25, 15, 100)
GetNetSolar <- function(Ra, Elev, th, tl, p){
  Vp = GetVp(p,th,tl)
  Rso <- (0.75+2*10^-5*Elev)*Ra
  Rs <- GetSolar(Ra, Elev, th, tl, p)
  Rnl <- 4.901*10^-9 * (1.35*Rs/(Rso+0.000001)-0.35) * (0.34 - 0.14 * Vp^0.5) * ((th+273.16)^4 + (tl+273.16)^4)/2
  Rns <- (1-0.23)*Rs
  Rn <- pmax(0,Rns - Rnl)
  return(Rn)}


#' Estimate plant transpiration component of potential evapotranspiration mediated by likely safe growing season temperatures
#'
#' @param th Mean daily high temperature (degrees Celsius)
#' @param tl Mean daily low temperature (degrees Celsius)
#'
#' @return Proportion of full potential evapotranspiration. Basically reduces the base formula for potential evapotranspiration to reflect plants may shut down their transpiration due to cold, leaving only evaporation from soil surfaces. This amount here is intended to be generalization for a diverse range of land cover, and is assumed 80% transpiration during growing season, and reduced to 0% during freezing conditions. Its actual amount is difficult to know, and would vary depending on plant cover.
#' @export
#'
#' @examples GetTransGrow(25, 15)
#' @examples GetTransGrow(10, -5)
#' @examples GetTransGrow(1, -15)
GetTransGrow <- function(th, tl) {#Adjust to reduction in transpiration due to cold, with evaporation only outside growing season
  ts = 0.8 #assumed T/ET ratio during growing season
  tw = 0 #assumed T/ET ratio during freezing season
  t <- (th+tl)/2
  tr <- 10 #generally as mean temperatures get below 10 transpiration shuts down, regardless of warm daytime temperatures
  G0 <- (t-0)/(tr)
  G1 <- pmin(1,pmax(0,G0)) #generally as mean temperatures get below 5 transpiration shuts down, regardless of warm daytime temperatures
  evmin = (tw)+(1-ts)
  G = G1*(1-evmin)+evmin
  return(G)}


GetDcl <- function(DayNumber){0.409*sin(2*3.141592*DayNumber/365-1.39)}


GetPET <- function(Ra, th, tl, p){
  Vpmax = 0.6108*exp(17.27*th/(th+237.3)) #saturation vapor pressure kPa
  Vpmin = 0.6108*exp(17.27*tl/(tl+237.3)) #saturation vapor pressure kPa
  logp <- log(p+1)
  e0 <- Ra*0.0508780  +
    Vpmax*0.7893714  +
    Vpmin*-0.5589255  +
    logp*-0.1309403  +
    Ra*Vpmax*0.0049383
  e <- pmax(0,e0)
  return(e)}


GetPETgs <- function(Ra, th, tl, p){
  e =  0.85*GetTransGrow(th, tl)*GetPET(Ra, th, tl, p) #0.85829 is crop coefficient to make comparable to Thornthwaite. Multiplied by growing season to zero out transpiration portion of evapotranspiration in freezing weather. Per day rate need to be multiplied by days per month to get monthly rate.
  return(e)}
