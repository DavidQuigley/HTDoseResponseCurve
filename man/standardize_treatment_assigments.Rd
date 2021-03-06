% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/synergy.R
\name{standardize_treatment_assigments}
\alias{standardize_treatment_assigments}
\title{Make two treatments consistently assigned to either treatment or treatment_2}
\usage{
standardize_treatment_assigments(D, treatment_1, treatment_2)
}
\arguments{
\item{D}{dataset}

\item{treatment_1}{String corresponding to the treatment to ensure is in 
the "treatment" column}

\item{treatment_2}{String corresponding to the treatment to ensure is in 
the "treatment_2" column}
}
\value{
A dataset with the same values as the dataset passed in as a 
parameter, but the treatment column will always be vehicle for treatment 1 or 
the value passed for treatment_1, and the treatment_2 column will always be 
the vehicle for treatment 2 or the value passed for treatment_2.
}
\description{
Will raise an error if there are more than two unique treatments in the 
dataset. Applied across all hours and plate IDs.
}

