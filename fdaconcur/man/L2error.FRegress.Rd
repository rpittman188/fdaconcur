\name{L2Error.fRegress}
\alias{L2Error.fRegress}
\title{L2Error.fRegress Function}
\usage{
Calculates the L2 difference across the curve using the given predictor and response datasets to predict the response curve and then compares the prediction to the true response curve that is known.
}
\description{
A L2 difference across each index using trapezoidal approximation. This function is typically not used by itself, rather, it is used inside the L2bestEst function to find the optimal number of basis functions to use when using the FRegress function.}
\examples{
NA #To be used with the L2bestEst function.
}