\name{InfluentialPoints}
\alias{InfluentialPoints}
\title{InfluentialPoints Function}
\usage{
InfluentialPoints(PredictorMat, ResponseMat, PredictorVec, BasisSelection = "L2", basislist, lambda = 10^(-1), plot = TRUE)
}
\arguments{
	PredictorMat is the matrix of predictor (X) values
	
	ResponseMat is the complete matrix of response (Y) values
	
	PredictorVec is the vector of X values that correspond to the missing, needed, Y values
	
	BasisSelection is the criteria used to select the number of Basis functions to use. ("L1" or "L2") which correspond to the "L1" and "L2" difference used in cross validation
	
	basisList is the user-derfined list of basis numbers that you want to optimize over. Must be odd for Fourier BasisType
	
	BasisType is the basis type. The two options are Fourier and BSpline.
	
	lambda is the smoothing parameter that can be found using the choosing_lambda function
	
	plot will give you several graphs and plots based on the influential points}

\description{
This function will show the user which events are the most influential to the final prediction of the missing response event. It starts by getting the prediction using all of the events, and then compares the prediction when each known event is left out. Then it tracks how much the prediction with one left out changes compared to the original prediction. This gives and indication of how influential that event is to the final result

}
\examples{
NA

}