\name{CVconcurrentmodel}
\alias{CVconcurrentmodel}
\title{CVconcurrentmodel Function}
\usage{
Leave1Outconcurrent(PredictorMat, ResponseMat, PredictorVec, BasisSelection = "L2", basislist, lambda = 10^(-1), plot = TRUE)
}
\arguments{
	PredictorMat is the matrix of predictor (X) values
	
	ResponseMat is the complete matrix of response (Y) values
	
	PredictorVec is the (X) vector for the main event where the corresponding (Y) is unknown
	
	BasisSelection is either "L1" or "L2" and is how the number of basis functions are selected for each column left out
	
	basisList is the user-derfined list of basis numbers that you want to optimize over.
	
	lambda is the smoothing parameter that can be found using the choosing_lambda function}

\description{
	This function leaves each column out and then runs the PredictFRegressNormTest function and predicts the (Y) values for the column that is left out. It then returns the same number of columns as the predictorMat and ResponseMat where each column is the prediction when that specific column is not included. It helps see how the model does predicting the known responses. Note that the data are selected to resemble an event that is not included in this analysis, so the results should be judge accordingly.}
\examples{
NA

}