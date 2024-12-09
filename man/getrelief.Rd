% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getrelief.R
\name{getrelief}
\alias{getrelief}
\title{Get topographic relief}
\usage{
getrelief(
  dm,
  r1 = NA,
  r2 = NA,
  n = 0,
  s = 0.1,
  p = c("low", "medium", "high", "exact"),
  breaks = NA,
  reverse = FALSE
)
}
\arguments{
\item{dm}{Digital elevation model (terra raster) with projected coordinates.}

\item{r1}{Minimum focal radius for estimating relief (only focal radius when n-0). Optional if using breaks.}

\item{r2}{Maximum focal radius for estimating relief. Optional if using breaks.}

\item{n}{Number of segments between maximum and minimum radii to estimate relief. n=0 for a simple focal range analysis using r1.}

\item{s}{Slope for which relief is estimated (Default is 0.1 for a 10\% slope). Should be a ratio of vertical units to horizontal units).}

\item{p}{Precision from low to exact with higher levels of precision requiring more processing time.}

\item{breaks}{Optional vector of relief thresholds to specifically calculate a radius rather than to depend on interpolation among a regular interval between maximum and minimum radii. Will combine with r1 and r2 vertical equivalence if they are populated.}

\item{reverse}{If true, give priority to widest focal neighborhoods, otherwise relief has to meet slope criteria for all the smaller focal neigborhoods first.}
}
\value{
Topographic relief of either a fixed radius or for a fixed slope for a range of radii. Relief above and below focal radius breaks is simple, not based on fixed slope.
}
\description{
Get topographic relief
}
\examples{
library(terra)#'
#load package data for elevation
data("denali")
dem <- denali
dem <- toraster(dem)
plot(dem)
#Reproject to get metric distance units.
dem <- reproject(dem, rs=250)

#Established desired vertical increment of focal range analysis.
prebreaks = c(50, 100,150,200,300,500,1000,2000,3000,4000,5000,6000,7000,8000)
maxrange <- max(minmax(dem)[2],1000)#Establish maximum relative to the highest elevation in DEM. This is to reduce the need for radii wider than required.
breaks <- prebreaks[prebreaks <=maxrange]
#Relief at a 10\% slope. Using "low" precision for faster performance
rng10 <- getrelief(dem, r1=1000, r2=maxrange*10, s=0.1, n=1, p='low', breaks = breaks)
plot(rng10)
#Relief at a 25\% slope. Using "medium" precision for better quality.
rng25 <- getrelief(dem, r1=1000, r2=maxrange*10, s=0.25, n=1, p='medium', breaks = breaks)
plot(rng25)

#get relief with a 10\% slope using 10 intervals between focal radii and default low precision.
rng <- getrelief(dem, r1=1000, r2=100000, s=0.1, n=10)
#get relief with a 25\% slope using 10 intervals between focal radii and default low precision.
rng <- getrelief(dem, r1=1000, r2=100000, s=0.25, n=10)
plot(rng)
#get 10\% relief using specific break points of vertical relief and medium precision.
rng <- getrelief(dem, r1=1000, r2=100000, s=0.1, n=1, p='medium', breaks = c(300, 1000, 2500, 5000))
plot(rng)
#get relief with a fixed radius of 2000 meters.
rng <- getrelief(dem, r1=2000, s=0.1, n=0)
plot(rng)

}
