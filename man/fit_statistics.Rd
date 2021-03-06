% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HTDoseResponseCurve.R
\name{fit_statistics}
\alias{fit_statistics}
\title{Fit all observations and return summary information from fit_stats data frame}
\usage{
fit_statistics(D, fct)
}
\arguments{
\item{D}{dataset}

\item{fct}{Non-linear function to fit, e.g. drc::LL.3(). See summary.}
}
\value{
data frame with fit summary data (is_fitted, hour, treatment, 
sample_type, fit_stats, ANOVA_F_test, ANOVA_P_value)
}
\description{
Fit all observations and return summary information from fit_stats data frame
}
\section{Non-linear function for curve fitting}{


Curve fitting is performed by the \code{drm()} function in the \code{drc} 
library. To fit the curve, you need to select a non-linear function. To 
estimate the slope, upper asymptote, lower asymptote, and EC50, pass 
drc::LL.4(). To fix the lower asymptote at 1 and estimate the other 
parameters, pass drc::LL.3(). To fix the upper asympotote at 1 and the lower 
asymptote at 0, pass drc::LL.2. For a list of available functions, see 
\code{drc::getMeanFunctions()}.
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
fit_statistics( ds, fct=LL.3() )
}
\seealso{
\code{\link{fit_DRC}}
}

