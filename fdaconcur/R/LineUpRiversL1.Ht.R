#Note that the specific height column is called Conght
#This is why we have MainRiver$CongHt

LineUpRiversL1.Ht<-function(MainRiver, NewCurve){
  library(pracma) #For trapezoid rule
  lower.index<-1 #starting
  upper.index<-nrow(NewCurve) #starting
  ResultsMat<-matrix(NA, nrow = 1, ncol =3)
  colnames(ResultsMat)<-c("LowerIndex", "UpperIndex", "L1Diff")
  selected.trimed.Cong<-NewCurve[lower.index:upper.index,]
  samesize.NewCurve<-approx(x=1:nrow(selected.trimed.Cong),y=selected.trimed.Cong$CongHt, n=nrow(MainRiver))$y

  absdiff<-abs(samesize.NewCurve - MainRiver$CongHt)*MainRiver$CongHt^2 #multiplying by this gives more weight at the peak of the curve.
  #I need the area under the absdiff curve
  L1Diff<-trapz(x=1:length(absdiff), y = absdiff)
  ResultsMat[1,]<-c(lower.index,upper.index,L1Diff)

  for(i in 1:50000){

    new.lower.index<-lower.index+1
    selected.trimed.Cong<-NewCurve[new.lower.index:upper.index,]
    samesize.NewCurve<-approx(x=1:nrow(selected.trimed.Cong),y=selected.trimed.Cong$CongHt, n=nrow(MainRiver))$y
    absdiff<-abs(samesize.NewCurve - MainRiver$CongHt)*MainRiver$CongHt^2
    L1Diff.loweradj<-trapz(x=1:length(absdiff), y = absdiff)
    L1lower.results<-c(new.lower.index,upper.index,L1Diff.loweradj)

    new.upper.index<-upper.index-1
    selected.trimed.Cong<-NewCurve[lower.index:new.upper.index,]
    samesize.NewCurve<-approx(x=1:nrow(selected.trimed.Cong),y=selected.trimed.Cong$CongHt, n=nrow(MainRiver))$y
    absdiff<-abs(samesize.NewCurve - MainRiver$CongHt)*MainRiver$CongHt^2
    L1Diff.upperadj<-trapz(x=1:length(absdiff), y = absdiff)
    L1upper.results<-c(lower.index,new.upper.index,L1Diff.upperadj)
    ResultsMat<-rbind(ResultsMat,L1lower.results,L1upper.results)

    #Idea is if adjusting the lower gives smaller L1 difference than adjusting upper, then change the lower
    #Otherwise change the upper
    lower.index<-ifelse(L1Diff.loweradj<L1Diff.upperadj,new.lower.index,lower.index)
    upper.index<-ifelse(L1Diff.loweradj<L1Diff.upperadj,upper.index,new.upper.index)

    if(upper.index-lower.index<10) break #Could be a bigger number than 10. just says to stop if we are not taking much of the new curve.
  }
  out<-list(ResultsMat[ResultsMat[,3]==(min(ResultsMat[,3])),], ResultsMat)
  names(out)<-c("Best Choice", "Results Matrix")
  return(out)
}
