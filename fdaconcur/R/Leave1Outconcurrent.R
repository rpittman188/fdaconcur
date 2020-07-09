#The leave one out and predict method is below:
Leave1Outconcurrent<-function(PredictorMat, ResponseMat, PredictorVec, BasisSelection = "L1", basislist, lambda = 10^(-1), plot = FALSE){
  AllPredictedResponse<-matrix(NA, ncol = ncol(ResponseMat), nrow = nrow(ResponseMat))
  AllPredictedLower<-matrix(NA, ncol = ncol(ResponseMat), nrow = nrow(ResponseMat))
  AllPredictedUpper<-matrix(NA, ncol = ncol(ResponseMat), nrow = nrow(ResponseMat))
  
  for(i in 1:ncol(PredictorMat)){ #FindBestBasis using either L1 or L2
    if(BasisSelection == "L1"){
      basisNum<-L1bestEst(PredictorMat[,-i], ResponseMat[,-i], basislist)$`Best Basis`
    }
    if(BasisSelection == "L2"){
      basisNum<-L2bestEst(PredictorMat[,-i], ResponseMat[,-i], basislist)$`Best Basis`
    }
    #tells it how many Fourier basis functions to use eachtime
    basisNum<-as.integer(names(basisNum))
    PredictedResponse<-PredictFRegressNormTest(PredictorMat[,-i], ResponseMat[,-i],PredictorVec , lambda = lambda, nBasis = basisNum) #Gets predicted Response for the ith column left out
    PredictedResponseEstimate<-PredictedResponse$PredictedResponse
    PredictedResponseLower<-PredictedResponse$Lower
    PredictedResponseUpper<-PredictedResponse$Upper
    AllPredictedResponse[,i]<-PredictedResponseEstimate #stores the actual prediction if needed
    AllPredictedUpper[,i]<-PredictedResponseUpper #stores the actual prediction if needed
    AllPredictedLower[,i]<-PredictedResponseLower #stores the actual prediction if needed
    
    #plots the true cong and cedar and the predicted cedar with intervals.
    if(plot == TRUE){
    par(mfrow=c(1,1))
    plot(x=Oct2015congDT, y= PredictorVec, lty= 2,type = "l", ylim = c(0,20), main = i)
    lines(x=Oct2015congDT,AllPredictedResponse[,i], lwd = 2,col = "red")
    #lines(x=Oct2015congDT,L1LinEst9, col = "purple")
    lines(x=Oct2015congDT,AllPredictedUpper[,i], col = "green")
    lines(x=Oct2015congDT,AllPredictedLower[,i], col = "green")
    lines(x=oct2015cedDT[1:96], y=Oct2015ced$Height[1:96],lwd=3, col="blue", lty=2) #Height we know
    lines(x=oct2015cedDT[97:240], y=Oct2015ced$Height[97:240],lwd=3, col="blue", lty=2) #Height we know
  }
  }
  out<-list(AllPredictedResponse, AllPredictedLower, AllPredictedUpper)
  names(out)<-c("Predicted Responses Without i", "Lower Bounds", "Upper Bounds")
  return(out)
}