#L2bestEst

#Make bestEstL2
L2bestEst<- function(xMat, yMat, basisList, BasisType = "Fourier", lambda = 10^-1){
  L2diff<- matrix(NA, ncol=ncol(xMat), nrow=length(basisList)) #Row is each basis trial, i is the test column for that itteration
  rownames(L2diff)<-as.character(basisList)
  for(i in 1:ncol(xMat)){
    trainX<-xMat[,-i]
    trainY<-yMat[,-i]
    testX<-xMat[,i]
    testY<-yMat[,i]
    for(j in 1:length(basisList)){
      L2diff[j,i]<-L2Error.fRegress(trainX, trainY, predVec=testX, trueY = testY, nBasis = basisList[j], Basis = BasisType, lambda = lambda)
    }
  }
  #after going through both loops all the times, we have a matrix with errors for all the inputed basis after using all of the columns as the test data.
  #now lets take the average of each column. Smallest means that basis has a smaller average L2 distance when leaving 1 out and using the rest to predict.
  average.L2diff<- apply(L2diff,1,mean)
  BestBasis<-which.min(average.L2diff)
  smallest.L2diff<-average.L2diff[BestBasis]
  out<-list(BestBasis,smallest.L2diff)
  names(out)<-c("Best Basis", "Smallest L2diff")
  return(out)
}
