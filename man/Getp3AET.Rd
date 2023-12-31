% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/XtremLow.R
\name{Getp3AET}
\alias{Getp3AET}
\title{Get peak 3-month evapotranspiration}
\usage{
Getp3AET(p, e)
}
\arguments{
\item{p}{vector of 12 monthly precipitation (mm)}

\item{e}{vector of 12 monthly potential evapotranspiration (mm)}
}
\value{
p3AET (mm)
}
\description{
Function calculates p3AET, the peak consecutive 3-month sum of the lesser among monthly precipitation and evapotranspiration. Ignoring soil water storage, this is equivalent to the expected actual evapotranspiration of the most productive 3 months. This is a measure of the degree to which precipitation is associated with warm temperatures. A Mediterranean climate with the same temperatures and precipitation as a subtropical savanna will tend to have a lower value than the savanna. Higher value would tend to associate with a greater proportion of precipitation as thunderstorms.
}
\examples{
mon <- 1:12
t <- generateTemp(-5,21)
#Wet summer
p <- generatePpt(700/12, 160, 5, 6, 12)
e <- GetPET(mon,43,t+5,t-5,p)
Getp3AET(p,e)
#Wet winter
p <- generatePpt(700/12, 160, 5, 12, 6)
e <- GetPET(mon,43,t+5,t-5,p)
Getp3AET(p,e)
}
