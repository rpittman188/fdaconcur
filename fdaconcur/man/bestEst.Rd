\name{bestEst}
\alias{bestEst}
\title{bestEst Function}
\usage{
bestEst(xMat, yMat, basisList, lambda, BasisType = "Fourier")
}
\arguments{
	xMat is the matrix of predictor (X) values
	
	yMat is the complete matrix of response (Y) values
	
	basisList is the user-derfined list of basis numbers that you want to optimize over.
	
	BasisType is the basis type. The two options are Fourier and BSpline.
	
	lambda is the smoothing parameter that can be found using the choosing_lambda function}

\description{
Gives the user the optimal number of basis functions to use when smoothing the data and calculating B0(t) and B1(t). Finds the MSE for every number of basis functions and types and returns the number that yields the lowest MSE based on the available data.}
\examples{
library(fda)
library(MASS)
library(corpcor)
hip1<-fda::gait
hip<-hip1[,, 'Hip Angle']
knee<-hip1[,, 'Knee Angle']
	#You could use BasisType = "BSpline" as well, and could use any odd number in the basisList.
outBestEst<-bestEst(hip[,-39], knee[, -39], lambda = 10^(-2), basisList = c(15,17,19,21), BasisType = "Fourier")

}