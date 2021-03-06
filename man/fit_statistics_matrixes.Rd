% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HTDoseResponseCurve.R
\name{fit_statistics_matrixes}
\alias{fit_statistics_matrixes}
\title{Convert fit statistics data frame to matrixes for further analysis}
\usage{
fit_statistics_matrixes(fits)
}
\arguments{
\item{fits}{data frame produced by call to \code{\link{fit_statistics}}}
}
\value{
list of lists, indexed first by sample type and then by 
ANOVA_P_value, AUC, and EC50. These elements are identically sized matrixes 
with the dimensions (treatment x hour)
}
\description{
This function is convenient for producing heat maps or other plots of 
fit summary data across a time series.
}
\examples{
sample_types = rep( c(rep("line1",3), rep("line2",3)), 5)
treatments = c(rep("DMSO",6), rep("drug",24))
concentrations = c( rep(0,6),rep(200,6), rep(500,6),rep(1000,6),rep(5000,6))
values=c(100,99,100,90,91,92,99,97,99,89,87,88,86,89,88,56,59,58,66,65,67,
         25,23,24,42,43,46,4,5,9)
hours = rep(48, length(values))
plate_id = "plate_1"
ds = create_dataset( sample_types=sample_types, treatments=treatments, 
                     concentrations=concentrations, hours=hours, 
                     values=values, plate_id=plate_id, 
                     negative_control = "DMSO")
library(drc)
fits = fit_statistics( ds, fct=LL.3() )
fit_statistics_matrixes( fits )
}

