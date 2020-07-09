\name{FindCedar}
\alias{FindCedar}
\title{FindCedar Function}
\usage{
FindCedar(StartDateTime, EndDateTime)
}
\arguments{
StartDateTime is the beginning time

EndDateTime is the end time}

\description{
	Used for the flood analysis to select the desired starting and stopping points of each event for the response matrix	
}
\examples{
	Need to run this first:
	AllCedar<-EntireCedarHt
AllCedar<-cbind(AllCedar,paste(AllCedar$Index, sep= " ", AllCedar$datetime))
names(AllCedar)<-c("Date", "Time", "Value", "Timestamp..UTC.05.00.")

FindCedar("1/1/95 1:00", "1/10/95 1:00"))

}