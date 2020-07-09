InfluentialPoints<-function(PredictorMat, ResponseMat, PredictorVec, BasisSelection = "L2", basislist, lambda = 10^(-1), plot = TRUE){
  #Decides how many basis functions to use for the full data
  if(BasisSelection == "L1"){
    basisNum<-L1bestEst(PredictorMat, ResponseMat, basislist)$`Best Basis`
  }
  if(BasisSelection == "L2"){
    basisNum<-L2bestEst(PredictorMat, ResponseMat, basislist)$`Best Basis`
  }
  basisNum<-as.integer(names(basisNum))
  FullPrediction<-PredictFRegress(PredictorMat,ResponseMat,PredictorVec, nBasis = basisNum, lambda = lambda)$PredictedResponse #gives the baseline with all included
  
  #gives the matrix of predictions 
  Leave1outPredictionMat<-Leave1Outconcurrent(PredictorMat,ResponseMat,PredictorVec, basislist = basislist, BasisSelection = BasisSelection)$`Predicted Responses Without i`
  
  #matrix of where the most influence comes from
  influence.Mat = FullPrediction - Leave1outPredictionMat
  
  influencescale.Mat<-scale(influence.Mat) #standardizes all of the differences
  
  #do some colors where darkred is < -2 darkgreen >2
  #red between -1 and -2 green between 1 and 2.
  #rest black
  #remember, positive means adding that event in will increase the prediciton. red means it decreases it
  
  colormat<-matrix(NA, nrow = nrow(influencescale.Mat), ncol = ncol(influencescale.Mat))
  for(i in 1:nrow(influencescale.Mat)){
    for(j in 1:ncol(influencescale.Mat)){
      if(influencescale.Mat[i,j]>2){
        colormat[i,j] = "Darkgreen"
      }
      else if(influencescale.Mat[i,j]<=2 & influencescale.Mat[i,j]>=1){
        colormat[i,j]="Green"
      }
      else if(influencescale.Mat[i,j]< -2 ){
        colormat[i,j]="Darkred"
      }
      else if(influencescale.Mat[i,j]>= -2 & influencescale.Mat[i,j]< -1){
        colormat[i,j]="Red"
      }
      else colormat[i,j]="Black"
    }
  }
  
  
  if(plot == TRUE){
  	par(mfrow=c(1,1))
    plot(influence.Mat[,1], ylim = c(-3,3),type = "l", ylab = "Height Difference", main = "Influence for Each Event")
    for(i in 2:ncol(influence.Mat)){
      lines(influence.Mat[,i], col = i)}
   #more plots 
    
    hist(influence.Mat, main = "Histogram of all Differences when 1 is left out")
    hist(influencescale.Mat, main ="Standardized Differences")
    
    #plots of the Each Column of the predictors with the color of how influential it is to the final prediction
    for(i in 1:ncol(influence.Mat)){
      plot(x=1:nrow(influence.Mat), y=PredictorMat[,i], cex = abs(influencescale.Mat[,i]), xlab = "Index", ylab = "Height", ylim = c(0,25),  col = colormat[,i], main = paste(i, "Event With Stardardized Influence"))
      legend( x="topright", 
              legend=c("Z < -2", "-2 <= Z < -1","-1 <= Z <= 1", "1 < Z <= 2", "Z > 2"), 
              col=c("darkred","red","black","green", "darkgreen"), lwd=1, 
              pch=c(16), merge=FALSE )
      }
  }
  out<-list(influence.Mat, influencescale.Mat, colormat, Leave1outPredictionMat, FullPrediction)
  names(out)<-c("Differences", "Standardized Differences", "ColorsMat", "Prediction When Left Out", "Full Prediction")
  return(out)
}