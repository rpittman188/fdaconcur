

library(MASS)
library(matrixcalc)
library(corpcor)
library(fda)
require(fda)

PredictFRegress <- function(PredictorMat,ResponseMat,predVec, nBasis, Basis="Fourier", lambda = 10^-1, Plot = FALSE){
  n<-nrow(PredictorMat)
  nPred<-ncol(PredictorMat) #Number of predictors (ie flood events)

  gaittime <- seq(1:n)
  gaitrange <- c(0,n)
  gaitfine = seq(0,n,.2)
  harmaccelLfd <- vec2Lfd(c(0, (2*pi/n)^2, 0), rangeval=gaitrange)
  if(Basis=="Fourier"){
    gaitbasis <- create.fourier.basis(gaitrange, nbasis=nBasis)} #original 15}

  if(Basis=="BSpline"){
    gaitbasis <- create.bspline.basis(gaitrange, nbasis=nBasis,norder=6) #original 20 leaving norder at 6.
  }

  mygaitExp <- array(NA, dim = c(n,nPred,2))
  mygaitExp[1:n, ,] <- seq(1:n)

  for(i in 1:nPred){
    mygaitExp[,i, 1] <- PredictorMat[,i]
    mygaitExp[,i, 2] <- ResponseMat[,i]
  }

  gaitSmooth = smooth.basisPar(gaittime, mygaitExp, gaitbasis, Lfdobj=harmaccelLfd, lambda=lambda) #Could Edit this part Original played around no different
  gaitfd = gaitSmooth$fd

  names(gaitfd$fdnames) = c("Normalized time", "Event", "Height")
  gaitfd$fdnames[[3]] = c("Cong", "Cedar")

  congfd  = gaitfd[,1] #Predictor Stuff
  cedfd = gaitfd[,2] # Response Stuff This seems like a more detailed version of what we get from the other textbook.

  xfdlist   = list(const=rep(1,nPred), cong=congfd) #number of observations here will change to 7
  betafdPar = fdPar(gaitbasis, harmaccelLfd)
  betalist  = list(const=betafdPar, cong=betafdPar)

  gaitRegress= fRegress(cedfd, xfdlist, betalist)

  # Figure 10.7
  op = par(mfrow=c(2,1))

  # Intercept
  betaestlist = gaitRegress$betaestlist
  cedIntercept = predict(betaestlist$const$fd, gaitfine)
  # dim(cedIntercept)

  # Cong coefficient
  congCoef = predict(betaestlist$cong$fd, gaitfine) #Slope Term

  cedhatfd = gaitRegress$yhatfd$fd
  cedhatmat = eval.fd(gaittime, cedhatfd)
  resmat. = mygaitExp[,,2] - cedhatmat #The 2 represents cedar
  SigmaE = cov(t(resmat.))

  cedfinemat   = eval.fd(gaitfine, cedfd)
  cedmeanvec   = eval.fd(gaitfine, mean(cedfd))
  cedhatfinemat= eval.fd(gaitfine, cedhatfd)
  resmat        = cedfinemat - cedhatfinemat
  ncurve        = dim(mygaitExp)[2] #Says the number of flood events
  resmat0 = cedfinemat - cedmeanvec %*% matrix(1,1,ncurve)
  SSE0 = apply((resmat0)^2, 1, sum)
  SSE1 = apply(resmat^2, 1, sum)
  ced.R2 = (SSE0-SSE1)/SSE0 #Truly R squared at each time point

  # Plot Cong Coefficient & Squared Multiple Correlation

  ylim2=c(0, max(congCoef, ced.R2)) #And this is what I had plotted before exactly.
  # Figure 10.8
  gaitbasismat = eval.basis(gaitfine, gaitbasis)
  y2cMap = gaitSmooth$y2cMap #Needed to get confidence intervals

  fRegressList1 = fRegress(cedfd, xfdlist, betalist,
                           y2cMap=y2cMap, SigmaE=SigmaE)

  fRegressList2 = fRegress.stderr(fRegressList1, y2cMap, SigmaE)
  betastderrlist = fRegressList2$betastderrlist

  error.B0<-predict(betastderrlist[[1]], gaitfine)
  error.B1<-predict(betastderrlist[[2]], gaitfine)

  if(Plot == TRUE){
    op = par(mfrow=c(2,1)) #These graphs do not work well in RStudio. Do in R.
    plotbeta(betaestlist, betastderrlist, seq(1,n,by=n/(n/10)))
    par(op)
    plot(gaitfine, congCoef, lwd=2, xlab='', ylab='', ylim=ylim2, type='l',
         main='Congaree Coefficient and Squared Multiple Correlation')
    abline(v=c(7.5, 14.7), lty='dashed')
    plot(gaitfine, ced.R2, type="l",main = "R-Squared") #Is this the R squared graph along with B1
    abline(h=1,col="Gray")
    plot(cedIntercept,type ="l",main = "Intercept Term")
    plot(congCoef,type="l",main="Slope Term")
  }

  yHat.t.mat<-matrix(NA,ncol = nPred , nrow = n)
  for(i in 1:nPred){
    yHat.t.mat[,i]<-approx(x=as.numeric(fRegressList1$yhatfdobj$argvals),y=fRegressList1$yhatfdobj$y[,i], n = n)$y
  }
  y.t<-ResponseMat

  MSE.t<-rep(NA,n)
  for(t in 1:n){
    SqError<-rep(NA,nPred)
    for(i in 1:nPred){
      SqError[i]<-(y.t[t,i]-yHat.t.mat[t,i])^2
    }
    MSE.t[t]<-sum(SqError)/(nPred-2)
  }

  epsilon<-matrix(NA, nrow = n, ncol = 1000)
  for(t in 1:n){
    epsilon[t,]<-rnorm(n=1000, mean = 0, sd = sqrt(MSE.t[t])) #1000 values for each t to add to function above
  }

  Xbar.t.predict<- predict(mean(congfd),gaitfine)
  samesize.Xbar<-approx(x=1:length(Xbar.t.predict), y=Xbar.t.predict, n = n)$y

  samesize.var.B0.t<-(approx(x=1:length(error.B0), y=error.B0, n = n)$y)^2
  samesize.var.B1.t<-(approx(x=1:length(error.B1), y=error.B1, n = n)$y)^2

  covB0B1.t<- rep(NA,n)
  for(t in 1:n){
    covB0B1.t[t]= -samesize.Xbar[t]*samesize.var.B1.t[t] #Needed since B0 and B1 are not independent
  }
  #Now we have a covariance of B0hat and B1hat for each t
  #We also have variance of B0 and variance of B1
  #Thus, we can create a variance covariance matrix at each time point
  #Will be 2x2x1220 where 1220 is n
  var.cov<- array(NA, dim = c(2,2,n), dimnames = list(c("B0","B1"),c("B0","B1"),rep("Time",n)))
  #So for each t 1:1220, we have var and cov of the B0hat and B1hat
  var.cov[1,1,]<-samesize.var.B0.t
  var.cov[2,2,]<-samesize.var.B1.t
  var.cov[1,2,]=var.cov[2,1,]=covB0B1.t
  #Now var.cov is an array with the variance covariance matrix at each time point t.
  #We also have cedIntercept congCoef

  #Before I actually initialize these, lets set them to be the same size as the rest.
  B0.t<-approx(x=1:length(cedIntercept), y=cedIntercept, n = n)$y
  B1.t<-approx(x=1:length(congCoef), y=congCoef , n = n)$y
  Slope = B1.t #Given in output
  Intercept = B0.t #Given in output

  predVec.samesize=approx(x=1:length(predVec), y=predVec, n = n)$y
  PredictedResponse = Intercept + predVec.samesize*Slope #Given in output

  #Now I can generate 1000 mean vectors at each time point. Then add my error at each time point to each
  #Need 1000 numbers at each time point 1:1220
  #Need to use whatever cong height I have.
  #This will be my X, in the function it will already be better looking.

  #For function: need the predVec to be same size as matrix
  #predVec.samesize

  generated.predictions<-matrix(NA, nrow = n, ncol = 1000) #each row is a time point t
  for(t in 1:n){ #goes through all time points good 1:750, 912-1220 is good
    b0b1Est<-matrix(0,nrow =1, ncol=2)
    b0b1Est<-mvrnorm(n=1000, mu = c(B0.t[t],B1.t[t]), Sigma = make.positive.definite(var.cov[,,t])) #Generates 1000 b0 and 1000 b1 at time t
    generated.predictions[t,] <- b0b1Est[,1]+b0b1Est[,2]*predVec.samesize[t]
  }
  generated.predictions.ep<-generated.predictions+epsilon
  #Now I have 1000 predictions for y at each t.

  upper.PI<-rep(NA,nrow(PredictorMat))
  lower.PI<-rep(NA,nrow(PredictorMat))
  for(t in 1:nrow(PredictorMat)){
    upper.PI[t]<-quantile(generated.predictions.ep[t,],.975)
    lower.PI[t]<-quantile(generated.predictions.ep[t,],.025)
  }

  out<-list(Intercept,Slope,PredictedResponse,upper.PI,lower.PI)
  names(out)<-c("Intercept","Slope","PredictedResponse","Upper","Lower")
  return(out)
}
