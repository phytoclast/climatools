% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dem_tools.R
\name{focalmax}
\alias{focalmax}
\title{Get focal maximum with specified radius}
\usage{
focalmax(x, r, p = c("low", "medium", "high", "exact"))
}
\arguments{
\item{x}{raster}

\item{r}{radius}

\item{p}{precision, from low to exact, with higher levels of precision requiring more processing time.}
}
\value{
focal maximum raster
This function makes use of lower resolution aggregated rasters to speed up calculations at increasing radii, compromising accuracy. Higher precision generates a focal neighborhood at a higher resolution, while lower precision aggregates more aggressively before focal analysis.
}
\description{
Get focal maximum with specified radius
}
\examples{
data("denali")
dem <- denali
dem <- toraster(dem)
dem <- reproject(dem, rs=250)
plot(dem)
x <- climatools::focalmax(dem, 2000)
plot(x)
}
