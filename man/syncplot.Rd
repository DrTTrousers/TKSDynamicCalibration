% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/output_functions.R
\name{syncplot}
\alias{syncplot}
\title{Syncronisation Plot}
\usage{
syncplot(x = tek_sum, title = "Peak Synchronisation Check")
}
\arguments{
\item{x}{TekScan summary from indexer() or autofit()}

\item{title}{A string to label the plot}
}
\value{

}
\description{
Graphically and interactively illustrates the synchronisation of TekScan output with Simulator output.
Useful graphical tool to check autofit() has correctly optimised the variables.
Best combined with modplot() as two halves of the same analysis, using plotly.
Can be partial bundled.
}
\examples{
\dontrun{
syncplot()
subplot(syncplot(), modplot(), nrows=2) \%>\% partial_bundle()
}
}
