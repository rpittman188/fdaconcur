
choosing_lambda<-function(dataX, dataY, nFourrierBasis = 15){
  require(fda)
n=nrow(dataX)

gaittime <- seq(1:n)
gaitrange <- c(0,n)
gaitfine = seq(0,n,1) #Changed from .2
harmaccelLfd1220 <- vec2Lfd(c(0, (2*pi/n)^2, 0), rangeval=gaitrange)
gaitbasis <- create.fourier.basis(gaitrange, nbasis=nFourrierBasis) #original 15

mygaitExp <- array(NA, dim = c(n,ncol(dataX),2))
mygaitExp[1:n, ,] <- seq(1:n)
for(i in 1:ncol(dataX)){
mygaitExp[,i, 1] <- dataX[,i]
mygaitExp[,i, 2] <- dataY[,i]
}

#Begin here assuming the mygaitExp is complete.
#####This part helps choose a negative lambda Shows that really anything below 6 is fine
gaitLoglam = seq(-10,25,0.25) #Maybe leave this.
nglam   = length(gaitLoglam)
gaitSmoothStats = array(NA, dim=c(nglam, 3),dimnames=list(gaitLoglam, c("log10.lambda", "df", "gcv") ) )
gaitSmoothStats[, 1] = gaitLoglam
for (ilam in 1:nglam) {
  gaitSmooth = smooth.basisPar(gaittime, mygaitExp, gaitbasis,
                               Lfdobj=harmaccelLfd1220, lambda=10^gaitLoglam[ilam])
  gaitSmoothStats[ilam, "df"]  = gaitSmooth$df
  gaitSmoothStats[ilam, "gcv"] = sum(gaitSmooth$gcv)}
# note: gcv is a matrix in this case
gaitSmoothStats
plot(gaitSmoothStats[, c(1, 3)], type='b') #This is GCV

#  set up plotting arrangements for one and two panel displays
#  allowing for larger fonts

op = par(mfrow=c(2,1))
par(op)
plot(gaitSmoothStats[, c(1, 3)], type="b", log="y") #so is this
plot(gaitSmoothStats[, 1:2], type="b", log="y") #THis is DF but intercept

}


