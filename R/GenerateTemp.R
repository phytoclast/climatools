#this generates a 12 month temperature curve from January and July temperatures based on simple sine wave.
#' Generate mean temperature monthly time series
#'
#' @param t01 January temperature (degrees Celsius)
#' @param t07 July temperature (degrees Celsius)
#'
#' @return Monthly time series of mean temperatures for 12 month period.
#' This function is useful to explore artificial climatic scenarios. Temperatures follow a simple sine wave with January and July as the coldest and warmest months (or the inverse for the southern hemisphere). A real climate may show some slight asymmetry between warm months and cold months, and be shifted in phase by a month or more, but this curve is very close to the most frequent pattern.
#' @export
#'
#' @examples generateTemp(-5,21)
generateTemp <- function(t01,t07){
  mon=c(1,2,3,4,5,6,7,8,9,10,11,12)
  t=round((t01-t07)/2*cos((mon-1)/12*2*3.141592)+(t01+t07)/2,1)
  return(t)}

#this generates a generalized precipitation curve from a mean, maximum, minimum,  precipitation, and maximum and minimum months.
#' Generate precipitation monthly time series
#'
#' @param p.mean mean monthly precipitation
#' @param p.max highest monthly precipitation
#' @param p.min lowest monthly precipitation
#' @param m.max number of month with highest precipitation
#' @param m.min number of month with lowest precipitation
#'
#' @return Monthly time series of precipitation for 12 month period.
#' This function is useful to explore artificial climatic scenarios. Precipitation follows a simple sine wave pattern if highest and lowest months are 6 months apart, and the mean precipitation is halfway between the highest and lowest values. Patterns other than a sine way will occur when high and low values are extreme or deviate from being 6 months apart. Only targeted mean monthly precipitation remains true in the latter case.
#' @export
#'
#' @examples generatePpt(p.mean=100, p.max=120, p.min=80, m.max=8, m.min=2)
generatePpt <- function(p.mean, p.max, p.min, m.max=7, m.min=1){
  mon=c(1,2,3,4,5,6,7,8,9,10,11,12)
  shift1 <- mon-m.min
  shift2 <- mon-m.max
  pshift1 <- (p.max-p.min)/2*cos((shift1+6)/12*2*3.141592)+(p.max+p.min)/2
  pshift2 <- (p.max-p.min)/2*cos((shift2)/12*2*3.141592)+(p.max+p.min)/2
  pwt1 <- ((cos((shift1)/12*2*3.141592)+1)/2)
  pwt2 <- ((cos((shift2)/12*2*3.141592)+1)/2)
  wts1 <- pwt1*(1-pwt2)
  wts2 <- pwt2*(1-pwt1)
  p0 <- ((wts1*pshift1)+(wts2*pshift2))/(wts1+wts2)
  exx<-seq(0.01,5,0.01)
  for(i in 1:length(exx)){#i=1
    dif0 = (mean(p0^exx[i]*max(p0)/max(p0^exx[i])) - p.mean)^2
    if(i==1){
      exdif = exx[1]
      dif=dif0
    }else{
      if(dif>=dif0){
        dif=dif0
        exdif = exx[i]
      }}}
  # rnd <- runif(12, 0.8,1.2)#option to add randomness to rain
  # p0 <- rnd*p0
  p1 <- p0^exdif*max(p0)/max(p0^exdif)
  p <- round(p1*p.mean/mean(p1),1)
  return(p)}
