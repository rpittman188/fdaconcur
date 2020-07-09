#Used for the flood analysis to select the desired starting and stopping points of each event for the response matrix

FindCedar<-function(StartDateTime, EndDateTime){
  begin.ind<-which(AllCedar$Timestamp..UTC.05.00. == StartDateTime)
  end.ind <-  which(AllCedar$Timestamp..UTC.05.00. == EndDateTime)
  range = end.ind - begin.ind #to make sure there are no missing number in the range
  CedarHt<-approx(x=1:nrow(AllCedar[begin.ind:end.ind,]), y=AllCedar[begin.ind:end.ind,"Value"], n = nrow(NewOct15Cong))$y
  dates<-AllCedar[begin.ind:end.ind,"Timestamp..UTC.05.00."]
  output<-list(CedarHt,dates,range)
  names(output)<-c("CedarHt", "Dates Returned" ,"Range")
  return(output)
}