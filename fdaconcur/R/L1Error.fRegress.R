#L1error.FRegress
#Will be used in bestEstL1.R
#Goal to calc L1 difference between true curve and predicted curve

  library(pracma)

L1Error.fRegress<-function(PredictorMat,ResponseMat,predVecX, trueY, nBasis, lambda = 10^-1, Basis="Fourier"){
  Predicted<-PredictFRegressNormTest(PredictorMat,ResponseMat,predVecX, nBasis = nBasis,Basis=Basis, lambda=lambda)$PredictedResponse
  True<-approx(x=1:length(trueY), y=trueY, n = length(Predicted))$y #Makes sure they are the same size and line up
  absdiff<-abs(True-Predicted)
  L1Diff<-trapz(x=1:length(absdiff), y = absdiff) #need library(pracma) 
  return(L1Diff) #L1 Diff for this particular basis number and specific singular test event
}
