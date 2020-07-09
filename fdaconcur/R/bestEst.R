bestEst<- function(xMat, yMat, basisList, BasisType = "Fourier", lambda = 10^-1){
  error<- matrix(NA, ncol=ncol(xMat), nrow=length(basisList)) #Row is each basis trial, i is the test column for that itteration
  rownames(error)<-as.character(basisList)
  for(i in 1:ncol(xMat)){
    trainX<-xMat[,-i]
    trainY<-yMat[,-i]
    testX<-xMat[,i]
    testY<-yMat[,i]
    for(j in 1:length(basisList)){
      error[j,i]<-Error.fRegress(trainX, trainY, predVec=testX, trueY = testY, nBasis = basisList[j], Basis = BasisType, lambda = lambda)
    }
  }
  #after going through both loops all the times, we have a matrix with errors for all the inputed basis after using all of the columns as the test data.
  #now lets take the average of each column
  average.error<- apply(error,1,mean)
  BestBasis<-which.min(average.error)
  smallest.error<-average.error[BestBasis]
  out<-list(BestBasis,smallest.error)
  names(out)<-c("Best Basis", "Smallest Error")
  return(out)
}
