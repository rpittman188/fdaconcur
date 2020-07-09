CVconcurrentmodel<-function(PredictorMat, ResponseMat, BasisSelection = "L2", basisList, lambda = 10^(-1), plot = TRUE)
{
  #Figures out which basis number to use each time
  PredictedResponse<-matrix(NA, ncol = ncol(ResponseMat), nrow = nrow(ResponseMat))

  for(i in 1:ncol(PredictorMat)){

  if(BasisSelection == "L1"){
    basisNum<-L1bestEst(PredictorMat[,-i], ResponseMat[,-i], basislist)$`Best Basis`
  }
  if(BasisSelection == "L2"){
    basisNum<-L2bestEst(PredictorMat[,-i], ResponseMat[,-i], basislist)$`Best Basis`
  }
  #tells it how many Fourier basis functions to use eachtime
  basisNum<-as.integer(names(basisNum))
  PredictedResponse[,i]<-PredictFRegressNormTest(PredictorMat[,-i], ResponseMat[,-i],PredictorMat[,i] , lambda = lambda, nBasis = basisNum)$PredictedResponse #Gets predicted Response for the ith column left out
  if(plot == TRUE){
  	par(mfrow=c(1,1))
    plot(ResponseMat[,i], type = "l", ylim = c(0,20),lwd = 2, ylab = "River Height (ft)", main = i)
    lines(PredictedResponse[,i], col = "red")
    legend( x="topright",
            legend=c("True Cedar Height", "Predicted Cedar Height"),
            col=c("black","red"), lwd=1)
    }
  }
  return(PredictedResponse)
}