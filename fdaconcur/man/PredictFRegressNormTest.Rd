\name{PredictFRegressNormTest}
\alias{PredictFRegressNormTest}
\title{Functional Data Analysis Using the Concurrent Model}
\usage{
PredictFRegressNormTest(PredictorMat,ResponseMat,predVec, nBasis, Basis="Fourier", lambda = 10^-1, Plot = FALSE, NormalErrors = FALSE) }
\arguments{
	PredictorMat is the matrix of predictor (X) values
	
	ResponseMat is the complete matrix of response (Y) values
	
	predVec is a vector of X values that correspond to the missing responce vector of interest.
	
	nBasis is the number of basis function needed to smooth the data. The optimal number can be found using the bestEst function from this package.
	
	Basis is the basis type. The two options are Fourier and BSpline.
	
	lambda is the smoothing parameter that can be found using the choosing_lambda function
	
	Plot gives you the option to produce plots for B0(t) and B1(t) along with their respective confidence intervals
	
	NormalErrors will return the p-values from the Shapiro-Wilks test of normality of errors at each index

}
\description{
This function produces an estimate for a missing vector y(t) given a known corresponding vector x(t). The function uses the known data to estimate B0(t) and B1(t) for the equation Y(t) = B0(t) + X(t)*B1(t), where Y(t) = ResponseMat and X(t) = PredictorMat. This function also produces a 95 percent prediction interval.}
\examples{
Using hip and knee data from the fda package.
library(fda)
library(MASS)
library(corpcor)
hip1<-fda::gait
hip<-hip1[,, 'Hip Angle']
knee<-hip1[,, 'Knee Angle']

pred39Knee<-PredictFRegress(hip[,-39], knee[, -39], predVec = hip[,39], lambda = 10^(-2), nBasis = 17 )

par( mfrow = c(1,1) )
plot(x= rownames(hip1[, , 'Knee Angle']), y = knee[, 39], type = "l", main = "Knee Angle", xlab = "Knee Angle", ylab = "Time")
lines(x= rownames(hip1[, , 'Knee Angle']), y = pred39Knee$PredictedResponse, col = "red" )
lines(x= rownames(hip1[, , 'Knee Angle']), y = pred39Knee$Lower, col = "green" )
lines(x= rownames(hip1[, , 'Knee Angle']), y = pred39Knee$Upper, col = "green" )
legend("topleft", legend=c("True Knee Angle", "Predicted Knee Angle", "95 PI"), col=c("black", "red", "green"), lty=1, cex=0.8)
}