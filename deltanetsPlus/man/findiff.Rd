% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/findiff.R
\name{findiff}
\alias{findiff}
\title{2nd order accuracy finite difference from three points}
\usage{
findiff(tpoints, f, tsi)
}
\arguments{
\item{tpoints}{Three consecutive time points in increasing order in an experiment}

\item{f}{Three columns of the log2FC matrix. Columns correspond to the samples given for \code{tpoints}.}

\item{tsi}{The indication of the current time point (1 for the initial time point, 2 for middle time points, and 3 for the last time point)}
}
\value{
\item{dfdt}{A vector of slopes at the current time point}
}
\description{
This function implements a 2nd order accurate finite difference for calculating slopes (time-derivatives) of the log2FC data using three time points. The function is used in \code{generateSlope} function.
}
\details{
For the initial time point (\code{tsi}=1) in an experiment, \code{tpoints} are the set of the initial and next two time ptoins. For the last time point in the experiment (\code{tsi}=3),\code{tpoints} are the set of the last and previous two time ptoins in increasing order. For a middle time point (\code{tsi}=2), \code{tpoints} are a current time point and two neighbor time points in increasing order.
Depending on \code{tsi}, three different strategies (forward, centeral, and backward) are used in approximation: forward for \code{tsi}=1, central for \code{tsi}=2, and backward for \code{tsi}=3.
}
