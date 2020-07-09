#Make a bestEstL1 Requires L1Error.fRegress but does it for every combination of basis number and event
L1bestEst<- function(xMat, yMat, basisList, BasisType = "Fourier",lambda=10^-1){
  L1diff<- matrix(NA, ncol=ncol(xMat), nrow=length(basisList)) #Row is each basis trial, i is the test column for that itteration
  rownames(L1diff)<-as.character(basisList)
  for(i in 1:ncol(xMat)){
    trainX<-xMat[,-i]
    trainY<-yMat[,-i]
    testX<-xMat[,i]
    testY<-yMat[,i]
    for(j in 1:length(basisList)){
      L1diff[j,i]<-L1Error.fRegress(trainX, trainY, predVec=testX, trueY = testY, nBasis = basisList[j], Basis = BasisType, lambda = lambda)
    }
  }
  #after going through both loops all the times, we have a matrix with errors for all the inputed basis after using all of the columns as the test data.
  #now lets take the average of each column. Smallest means that basis has a smaller average L1 distance when leaving 1 out and using the rest to predict.
  average.L1diff<- apply(L1diff,1,mean) #each row mean which is for each specific number of basis functions in the basisList vector
  BestBasis<-which.min(average.L1diff)
  smallest.L1diff<-average.L1diff[BestBasis]
  out<-list(BestBasis,smallest.L1diff)
  names(out)<-c("Best Basis", "Smallest L1diff")
  return(out)
}
