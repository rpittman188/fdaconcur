\name{choosing_lambda.fun}
\alias{choosing_lambda.fun}
\title{choosing_lambda.fun Function}
\usage{
This function uses the provided data and produces graphs that show the optimal smoothing parameter to use.

choosing_lambda(dataX, dataY)
}
\description{
This function uses the provided data and produces graphs that show the optimal smoothing parameter to use. The process uses generalized cross valdidation (gcv) at a predetermined grid of lambdas (-6,6). The optimal lambda is 10^(i) where i is the log10.lambda that produces the lowest gcv.}
\examples{
	#The duration of the example uses hip and knee data from the fda package. It uses 38 complete known hip/knee angle relationships to predict a 39th. The prediction can then be compared to the known knee angle which was not used in the calculation.
library(fda)
library(MASS)
library(corpcor)
hip1<-fda::gait
hip<-hip1[,, 'Hip Angle']
knee<-hip1[,, 'Knee Angle']

choosing_lambda(hip[,-39], knee[, -39])

}