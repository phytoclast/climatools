% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/XtremLow.R
\name{Biotemperature}
\alias{Biotemperature}
\title{Mean annual biotemperature}
\usage{
Biotemperature(t)
}
\arguments{
\item{t}{Vector of monthly mean temperatures for 12 months (degrees Celsius)}
}
\value{
Mean annual temperature between 0 and 30 (degrees Celsius; all monthly temperatures below zero counted as zero, and those above 30 limited to 30, prior to averaging)
}
\description{
Mean annual biotemperature
}
\examples{
t <- generateTemp(-5,21)
Biotemperature(t)
}
\references{
Holdridge, L.R. and Tosi Jr, J.A., 1967. Tropical Science Center: San jose, Costa Rica. Life Zone Ecology; Tropical Science Center: San Jose, Costa Rica.
}
