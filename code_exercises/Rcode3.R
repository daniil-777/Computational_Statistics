# R-code to illustrate confidence and prediction intervals 
#   and bias-variance trade-off. 

###########################################################################
###########################################################################

# Confidence intervals and prediction intervals:
install.packages("ISwR")

library(ISwR)  # Introductory Statistics with R (Dalgaard)
data(thuesen)
? thuesen
summary(thuesen)

# remove missing value:
thuesen <- thuesen[!is.na(thuesen[,2]),]
summary(thuesen)

# Create first plot:
plot(short.velocity ~ blood.glucose, data=thuesen)

# Fit simple linear regression
fit <- lm(short.velocity ~ blood.glucose, data=thuesen)
summary(fit)
abline(fit)

# 95% CI for intercept by hand:
n <- nrow(thuesen)
coef(fit)[1] - qt(.975, n-2)*.11748
coef(fit)[1] + qt(.975, n-2)*.11748

# Automatic:
confint(fit)

#############################################################

# Next, we want a confidence intervals for E(y_0), 
#   corresponding to some new value x0=10:
# We are interested in the *average* velocity for patients with 
#   a blood glucose of 10

x0 <- 10
(pred.frame <- data.frame(blood.glucose=x0))

# Obtain fitted value:
predict(fit, pred.frame)

# Check by hand:
(fitted <- fit$coef[1] + fit$coef[2]*x0)

# Obtain confidence interval for E(y0):
predict(fit, pred.frame, level=.95, interval="c")

# Check by hand:
# Quantile of t-distribution that we need:
(quant <- qt(.975,n-2))
# Sigma.hat value:
(sigma.hat <- sqrt(sum((fit$resid)^2)/(n-2)))
# Design matrix:
X <- as.matrix(cbind(1,thuesen[,1]))
XtXi <- solve(t(X) %*% X)
x00 <- as.matrix(c(1,x0),nrow=2)
# Now compute estimate for sd(\hat y0):
(se <- sigma.hat * sqrt( t(x00) %*% XtXi %*% x00))
# Confidence interval:
(lower <- fitted - quant*se)
(upper <- fitted + quant*se)

# Double check:
predict(fit, pred.frame, level=.95, interval="c")

###############################################################
  
# We now want a confidence interval for y_0, 
#   also called prediction interval. 
# We are interested in the velocity for a patient with 
#   a blood glucose of 10

# By hand:
(se <- sigma.hat * sqrt( 1 + t(x00) %*% XtXi %*% x00))
(lower <- fitted - quant*se)
(upper <- fitted + quant*se)

# Automatically:
predict(fit, pred.frame, level=.95, interval="p")

###########################################################

# Finally, let's consider the *pointwise* confidence/prediction band 

# Use x-values 4,5,...,20
x <- c(4:20)
(pred.frame <- data.frame(blood.glucose=x))
CI <- predict(fit, pred.frame, interval="c")
PI <- predict(fit, pred.frame, interval="p")

# Create plot. We use ylim to make sure that 
# the confidence&prediction intervals fit in the plot
plot(short.velocity ~ blood.glucose, data=thuesen, 
          ylim=range(PI,thuesen$short.velocity))
matlines(x, CI, lty=c(1,2,2), col="red", lwd=2)
matlines(x, PI, lty=c(1,3,3), col="blue", lwd=2)

# What is the behavior of the CIs as n grows?
# What is the behavior of the PIs as n grows?



###########################################################################
###########################################################################

# R code to illustrate bias-variance trade-off

# create non-linear function:
f <- function(x){
  .3* x - 0.2*x^2 + 0.1*x^3 + sin(2*x) 
}

# plot function:
grid <- seq(from=-5,to=5,length=300)
par(mfrow=c(1,1))
plot(grid,f(grid),type="l")

# we are going to use loess curves 
#   (loess=locally estimated scatterplot smoothing)
#   with different levels of smoothing:
span <- c(0.1,0.2,0.3,0.45,0.7,1)    # smoothing parameter for loess:
# small values mean little smoothing and hence complex models 

# parameters for data simulation:
sigma <- 1.5                        # standard devation of noise
n <- 100                            # sample size
x <- seq(from=-5,to=5, length=n)    # x-values (fixed throughout simulation)

# simulate data:
y <- f(x) + rnorm(n=length(x),mean=0,sd=sigma)

# fit loess smoothers and plot the data and the fitted regressions:
par(mfrow=c(2,3))
for (i in 1:length(span)){
  plot(x,f(x), type="l", lwd=2, main=paste("alpha=",span[i]))
  points(x,y)
  lo <- loess(y ~ x, span=span[i])
  lines(grid, predict(object=lo, grid),col="red")
}

# specify xtest, i.e., the x value of your test point.
#   (on the blackboard we called this x0)
xtest <- -2

# repeat the above 20 times. 
# (we now only plots the fitted regression, not the data sets.) 
par(mfrow=c(2,3))
for (i in 1:length(span)){
  plot(x,f(x), type="l", lwd=2, main=paste("alpha=",span[i]))
  for(j in 1:20){
    y <- f(x) + rnorm(n=length(x),mean=0,sd=sigma)
    lo <- loess(y ~ x, span=span[i])
    lines(x, predict(object=lo, x),col="gray")
    abline(v=xtest, lty=3)
  }
}

##################################################

# larger simulation

# parameters:
nsim <- 1000            # number of simulations

# create objects to store fitted value at xtest:
fit.test <- matrix(rep(NA,nsim*length(span)),nrow=nsim)
#   rows correspond to different simulations, 
#   columns correspond to different smoothing parameters

# create object to store true y value at xtest:
y.test <- rep(NA,nsim)

for (i in 1:nsim){
  # generate new training data:
  y <- f(x) + rnorm(n=length(x),mean=0,sd=sigma)
  # generate new test data and store it in ith element of y.test
  y.test[i] <- f(xtest) + rnorm(n=1,mean=0,sd=sigma)
  
  # fit loess curves with different smoothing levels
  for (j in 1:length(span)){ 
    lo <- loess(y ~ x, span=span[j])
    # store fitted value at xtest: 
    fit.test[i,j] <- predict(object=lo, xtest)
  }
} 

# plot histograms of fitted values at xtest:
par(mfrow=c(2,3))
for (i in 1:length(span)){
  hist(fit.test[,i],xlim=range(fit.test),freq=F,
       main=paste("alpha=",span[i]),xlab=paste("fitted value at x=",xtest))
  lines(density(fit.test[,i]))
  abline(v=f(xtest), col="red")
}

# Check bias variance decomposition:
(ExpTestMSE <- apply((fit.test-y.test)^2,2,mean))
(Bias2 <- (apply(fit.test,2,mean)-f(xtest))^2)
(Var <- apply(fit.test,2,var))
(VarEps <- var(y.test))
Bias2+Var+VarEps - ExpTestMSE
# Note small errors because cross-terms do not fully disappear in simulation

# Plot results:
par(mfrow=c(1,1))
plot(span,ExpTestMSE,ylim=range(0,max(ExpTestMSE)), 
     type="l",lwd=2,col="purple", xlab="smoothing parameter alpha (small values indicate a flexible model)", ylab="")
lines(span, Bias2,lwd=2,col="red")
lines(span,Var,lwd=2,col="blue")
abline(h=VarEps,lty=2,lwd=2)
#legend("topleft", c("Expected test MSE at xtest", "Bias^2", "Variance", "Irreducible error"),lwd=2, lty=c(1,1,1,2),col=c("purple", "red", "blue","black"), bty="n", cex=0.5)



