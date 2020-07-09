#L2error.FRegress
#Will be used in L2bestEst.R
#Goal to calc L2 difference between true curve and predicted curve

 library(pracma)
L2Error.fRegress<-function(PredictorMat,ResponseMat,predVecX, trueY, nBasis, lambda = 10^-1, Basis="Fourier"){
  Predicted<-PredictFRegressNormTest(PredictorMat,ResponseMat,predVecX, nBasis = nBasis,Basis=Basis, lambda=lambda)$PredictedResponse
  True<-approx(x=1:length(trueY), y=trueY, n = length(Predicted))$y #Makes sure they are the same size and line up
  squaredifference<-(True - Predicted)^2
  L2Diff<-sqrt(trapz(x=1:length(squaredifference), y = squaredifference)) #need library(pracma) 
  return(L2Diff) #L2 Difference/Distance
}
