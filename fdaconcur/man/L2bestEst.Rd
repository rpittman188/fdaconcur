\name{L2bestEst}
\alias{L2bestEst}
\title{L2bestEst Function}
\usage{
L2bestEst(xMat, yMat, basisList, lambda=10^-1, BasisType = "Fourier")
}
\arguments{
	xMat is the matrix of predictor (X) values
	
	yMat is the complete matrix of response (Y) values
	
	basisList is the user-derfined list of basis numbers that you want to optimize over. (If using Fourier it must be odd numbers)
	
	BasisType is the basis type. The two options are Fourier and BSpline.
	
	lambda is the smoothing parameter that can be found using the choosing_lambda function}

\description{
Gives the user the optimal number of basis functions to use when smoothing the data and calculating B0(t) and B1(t). Finds the average L2 distance for every number of basis functions and types then returns the number of basis that yields the lowest L2 distance between predicted curves and the true curves based on the available data.}
\examples{
library(fda)
library(MASS)
library(corpcor)
hip1<-fda::gait
hip<-hip1[,, 'Hip Angle']
knee<-hip1[,, 'Knee Angle']
	#You could use BasisType = "BSpline" as well, and could use any odd number in the basisList.
outBestEst<-L2bestEst(hip[,-39], knee[, -39], lambda = 10^(-2), basisList = c(15,17,19,21), BasisType = "Fourier")

}