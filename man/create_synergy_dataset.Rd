% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HTDoseResponseCurve.R
\name{create_synergy_dataset}
\alias{create_synergy_dataset}
\title{Create a synergy dataset from raw data without a plate map}
\usage{
create_synergy_dataset(sample_types, treatments_1, treatments_2,
  concentrations_1, concentrations_2, values, hours = 0,
  plate_id = "plate_1", negative_control = NA, summary_method = "mean")
}
\arguments{
\item{sample_types}{vector of sample types}

\item{treatments_1}{vector of treatment one}

\item{treatments_2}{vector of treatment two}

\item{concentrations_1}{vector of concentrations of treatment one}

\item{concentrations_2}{vector of concentrations of treatment two}

\item{values}{vector of measured response to treatments one and two}

\item{hours}{time points for each observation. If a number, the same time 
point is assigned to all observations. If a vector, there should be one 
number for each observation. Defaults to 0.}

\item{plate_id}{experiment identification string, useful if multiple datasets 
are later combined. Defaults to "plate_1"}

\item{negative_control}{A designation for the negative controls in this 
dataset, if they exist. Value may be NA, a number, a string, or a data frame.

\itemize{
 \item{NA: Use when there are no negative control measurements. The contents 
   of the column named 'value_normalized' will be copied from the contents of 
   the column named 'value'. }
 \item{Number: Use when each treatment has been labeled with a concentration 
   (typically 0) that indicates the vehicle control. Each treatment must 
   contain one or more observations with this concentration, and these 
   observations will be the negative controls.}
 \item{string: Use when a single set of observations is a universal control. 
   The treatment whose name matches the string is the universal 
   negative control all of the data.}
 \item{data frame: Use when more than one negative control exists, and you 
   have to map different treatments to a particular negative control. The 
   data frame must have names 'drug' and 'vehicle', and the data frame will 
   map match treatments in the 'drug' column to those in the 
   'vehicle' column.}
}}

\item{summary_method}{summarize replicate measures by either mean or median; 
must be one of "mean", "median". Defaults to "mean"}
}
\value{
a data frame where columns indicate the sample type, treatment 1, 
treatment 2, concentration of treatment 1, concentration of treatment 2, 
observed raw value, normalized value, name of the negative_control treatment, 
whether a particular row is a negative control for at least one other row, 
hours since the start time, and plate of origin.
}
\description{
A synergy dataset differs from a standard dataset in that each value is the 
result of combining two distinct treatments in two distinct concentrations. 
If there are measurements where only one treatment was present, the other 
treatment should be specified to have a concentration equal to zero.
}
\examples{
dose_SCH=c(0.1, 0.5, 1, 2, 4)
eff_SCH = c(0.6701, 0.6289, 0.5577, 0.4550, 0.3755)
dose_4HPR = c(0.1, 0.5, 1, 2)
eff_4HPR = c(0.7666, 0.5833, 0.5706, 0.4934)
eff_comb = c(0.6539, 0.4919, 0.3551, 0.2341)
syn = data.frame( 
    treatment_1 = rep("SCH66336", 13),
    conc_1 = c( dose_SCH, rep(0, 4), dose_SCH[1:4]),
    treatment_2 = rep("4-HPR", 13),
    conc_2 = c( rep(0, 5), dose_4HPR, dose_4HPR ),
    values = c(eff_SCH, eff_4HPR, eff_comb ), stringsAsFactors=FALSE )
ds_lk = create_synergy_dataset( sample_types = rep("sample_1", 13), 
                                treatments_1 = syn$treatment_1,
                                treatments_2 = syn$treatment_2,
                                concentrations_1 = syn$conc_1,
                                concentrations_2 = syn$conc_2,
                                values = syn$values)
}

