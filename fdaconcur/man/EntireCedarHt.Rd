\name{EntireCedarHt}
\alias{EntireCedarHt}
\title{EntireCedarHt Data}
\usage{
This is raw Cedar Creek Height data from 1/1/95 to 10/1/19. There are 671,818 records in mostly 15 minute intervals. Several intervals of time are missing, including October 2015.
}
\description{
The data is from:
Jason M. Fine | U.S. Geological Survey
Acting Associate Director for Data - Caribbean Florida Water Science Center 
Assistant Director for Data - North Carolina - South Atlantic Water Science Center 
On February 4, 2020.
}
\examples{
Run this to use the FindCedar Function:
AllCedar<-EntireCedarHt
AllCedar<-cbind(AllCedar,paste(AllCedar$Index, sep= " ", AllCedar$datetime))
names(AllCedar)<-c("Date", "Time", "Value", "Timestamp..UTC.05.00.")

}