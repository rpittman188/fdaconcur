\name{Leave1Outconcurrent}
\alias{Leave1Outconcurrent}
\title{Leave1Outconcurrent Function}
\usage{
Leave1Outconcurrent(PredictorMat, ResponseMat, PredictorVec, BasisSelection = "L1", basislist, lambda = 10^(-1), plot = FALSE)
}
\arguments{
	PredictorMat is the matrix of predictor (X) values
	
	ResponseMat is the complete matrix of response (Y) values
	
	PredictorVec is the (X) vector for the main event where the corresponding (Y) is unknown
	
	BasisSelection is either "L1" or "L2" and is how the number of basis functions are selected
	
	basisList is the user-derfined list of basis numbers that you want to optimize over.
	
	lambda is the smoothing parameter that can be found using the choosing_lambda function}

\description{
	This function leaves each column out and then runs the PredictFRegressNormTest function and predicts the (Y) values from the PredictorVec. It then returns the same number of columns as the predictorMat and ResponseMat where each column is the prediction when that specific column is not included. It helps see what happens to the final prediction when each event is not included to help see how influential it is to the final prediction.
}
\examples{
NA

}