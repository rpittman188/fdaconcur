LineUpRiversL2<-function(MainRiver, NewCurve){
  library(pracma) #For trapezoid rule
  lower.index<-1 #starting
  upper.index<-nrow(NewCurve) #starting
  ResultsMat<-matrix(NA, nrow = 1, ncol =3)
  colnames(ResultsMat)<-c("LowerIndex", "UpperIndex", "L2Diff")
  selected.trimed.Cong<-NewCurve[lower.index:upper.index,]
  samesize.NewCurve<-approx(x=1:nrow(selected.trimed.Cong),y=selected.trimed.Cong$CongHt, n=nrow(MainRiver))$y

  squaredifference<-(samesize.NewCurve - MainRiver$CongHt)^2
  #I need the area under the squaredifference curve
  L2Diff<-sqrt(trapz(x=1:length(squaredifference), y = squaredifference))
  ResultsMat[1,]<-c(lower.index,upper.index,L2Diff)

  for(i in 1:50000){

    new.lower.index<-lower.index+1
    selected.trimed.Cong<-NewCurve[new.lower.index:upper.index,]
    samesize.NewCurve<-approx(x=1:nrow(selected.trimed.Cong),y=selected.trimed.Cong$CongHt, n=nrow(MainRiver))$y
    squaredifference<-(samesize.NewCurve - MainRiver$CongHt)^2
    L2Diff.loweradj<-sqrt(trapz(x=1:length(squaredifference), y = squaredifference))

    L2lower.results<-c(new.lower.index,upper.index,L2Diff.loweradj)

    new.upper.index<-upper.index-1
    selected.trimed.Cong<-NewCurve[lower.index:new.upper.index,]
    samesize.NewCurve<-approx(x=1:nrow(selected.trimed.Cong),y=selected.trimed.Cong$CongHt, n=nrow(MainRiver))$y
    squaredifference<-(samesize.NewCurve - MainRiver$CongHt)^2
    L2Diff.upperadj<-sqrt(trapz(x=1:length(squaredifference), y = squaredifference))
    L2upper.results<-c(lower.index,new.upper.index,L2Diff.upperadj)

    ResultsMat<-rbind(ResultsMat,L2lower.results,L2upper.results)

    #Idea is if adjusting the lower gives smaller L2 difference than adjusting upper, then change the lower
    #Otherwise change the upper
    lower.index<-ifelse(L2Diff.loweradj<L2Diff.upperadj,new.lower.index,lower.index)
    upper.index<-ifelse(L2Diff.loweradj<L2Diff.upperadj,upper.index,new.upper.index)

    if(upper.index-lower.index<10) break #Could be a bigger number than 10. just says to stop if we are not taking much of the new curve.
  }
  out<-list(ResultsMat[ResultsMat[,3]==(min(ResultsMat[,3])),], ResultsMat)
  names(out)<-c("Best Choice", "Results Matrix")
  return(out)
}
