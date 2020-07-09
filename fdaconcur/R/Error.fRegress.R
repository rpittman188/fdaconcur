#My Error Function
Error.fRegress<-function(PredictorMat,ResponseMat,predVecX, trueY, nBasis, lambda = 10^-1, Basis="Fourier"){
  Predicted<-PredictFRegress(PredictorMat,ResponseMat,predVecX, nBasis = nBasis,Basis=Basis, lambda=lambda)$PredictedResponse
  True<-approx(x=1:length(trueY), y=trueY, n = length(Predicted)) #Makes sure they are the same size and line up
  ErrorVec<-rep(NA,length(Predicted))
  for(i in 1:length(Predicted)){
    ErrorVec[i]<- (Predicted[i]-True$y[i])^2 #length n
  }
  MSE<-mean(ErrorVec) #MSE
  return(MSE) #MSE
}
