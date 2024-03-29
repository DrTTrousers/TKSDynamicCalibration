% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/timepoint_functions.R
\name{data_join}
\alias{data_join}
\title{Data Joining}
\usage{
data_join(r = results, s = simdata)
}
\arguments{
\item{r}{Calibrated Tekscan results from sumdataframe(), or list.}

\item{s}{Simulator data as extracted by simext_full()}
}
\value{

}
\description{
A function to combine the results with the simdata for use plotting the simulator peaks.
Currently only used for simulator peaks, but could be useful for regression analysis.
Only combines data where the simulator and recording have the same names.
}
\examples{
\dontrun{
data_join()
}
}
