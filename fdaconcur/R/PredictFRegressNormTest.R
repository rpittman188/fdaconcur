#This function uses functional data analysis and the concurrent model to predict the level of a missing curve using the known Xt curve. It also gives a prediction interval based on normally distributed errors and checks that that assumption is valid.

library(MASS)
#library(matrixcalc) take out see if things can still run
library(corpcor)
library(fda)
require(fda)


PredictFRegressNormTest <- function(PredictorMat,ResponseMat,predVec, nBasis, Basis="Fourier", lambda = 10^-1, Plot = FALSE, NormalErrors = FALSE){
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
  
  xfdlist   = list(const=rep(1,nPred), cong=congfd) 
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
  
  cedhatfd = gaitRegress$yhatfd$fd #
  cedhatmat = eval.fd(gaittime, cedhatfd)  #1 to 2000
  resmat. = mygaitExp[,,2] - cedhatmat #The 2 represents cedar (response)
  SigmaE = cov(t(resmat.)) #Estimated after the prelim anyalysis with fRegress
  
  cedfinemat   = eval.fd(gaitfine, cedfd)
  cedmeanvec   = eval.fd(gaitfine, mean(cedfd))
  cedhatfinemat= eval.fd(gaitfine, cedhatfd)
  resmat        = cedfinemat - cedhatfinemat #finer grid than the .
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
  
  #Main analysis that uses SigmaE
  fRegressList1 = fRegress(cedfd, xfdlist, betalist,
                           y2cMap=y2cMap, SigmaE=SigmaE)
  
  fRegressList2 = fRegress.stderr(fRegressList1, y2cMap, SigmaE)
  betastderrlist = fRegressList2$betastderrlist
  
  error.B0<-predict(betastderrlist[[1]], gaitfine)
  error.B1<-predict(betastderrlist[[2]], gaitfine)
  
  if(Plot == TRUE){
    op = par(mfrow=c(1,1)) #These graphs do not work well in RStudio. Do in R.
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
  
  #Predict for all the values flood events
  yHat.t.mat<-matrix(NA,ncol = nPred , nrow = n)
  #takes the predicted heights for the other known events and interpolates to make them the same size (n)
  #Note these predictions use the same number of basis functions as we use for the main event
  #Other function in the package can make a more appropriate actual prediction for CV purposes
  for(i in 1:nPred){
    yHat.t.mat[,i]<-approx(x=as.numeric(fRegressList1$yhatfdobj$argvals),y=fRegressList1$yhatfdobj$y[,i], n = n)$y
  }
  y.t<-ResponseMat
  
  #test for normal residuals after the main analysis. Makes sure we are okay to use the normal distribution to get prediction interval
      updatedresmat. = mygaitExp[,,2] - yHat.t.mat #mygaitExp[,,2] same as y.t gets y-yHat at each index
  if(NormalErrors == TRUE){
    #Now want to check that these residuals are normal
    sw.pval<-c()
    for(i in 1:nrow(updatedresmat.)){
      sw.pval[i]<-shapiro.test(updatedresmat.[i,])$p.value #does test for normality at every point
    }
  }
  
  MSE.t<-rep(NA,n)
  for(t in 1:n){
    SqError<-rep(NA,nPred)
    for(i in 1:nPred){
      SqError[i]<-(updatedresmat.[t,i])^2 # at each time point for each flood was : (y.t[t,i]-yHat.t.mat[t,i])^2
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
  #Need 1000 numbers at each time point 1:2000
  #Need to use whatever cong height I have.
  #This will be my X, in the function it will already be better looking.
  
  #For function: need the predVec to be same size as matrix
  #predVec.samesize
  
  generated.predictions<-matrix(NA, nrow = n, ncol = 1000) #each row is a time point t
  for(t in 1:n){ 
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
  
  if(NormalErrors == TRUE){ #returns info about normal errors
    out<-list(Intercept,Slope,PredictedResponse,upper.PI,lower.PI, updatedresmat., sw.pval, yHat.t.mat)
    names(out)<-c("Intercept","Slope","PredictedResponse","Upper","Lower", "General Residual Matrix", "Shap-Wilks p-values", "yHat Matrix")
  }
  
  else{
  out<-list(Intercept,Slope,PredictedResponse,upper.PI,lower.PI, yHat.t.mat)
  names(out)<-c("Intercept","Slope","PredictedResponse","Upper","Lower", "yHat Matrix") }
  
  return(out)
}