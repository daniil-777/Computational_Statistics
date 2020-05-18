# drago.plecko@stat.math.ethz.ch

set.seed(21)
nrep <- 100
slope <- numeric(nrep)

x <- rnorm(40, 20, 3)
for (i in 1:nrep){
  ##y <- 1 + 2 * x + 5 * rnorm(length(x)) 
  y   <- 1 + 2 * x + 5 * rnorm(length(x), mean = 0, sd = (x-15)^2 / 30)
  reg <- lm(y ~ x)
  slope[i] <- coefficients(reg)[2]
}
  
X <- cbind(1, x)
XtXinv <- solve(crossprod(X)) ## crossprod(X) <=> t(X) %*% X = X'X > tvar <- 5^2 * XtXinv[2, 2]
tvar <- 5^2 * XtXinv[2, 2]
  
hist(slope, breaks = 15, freq = FALSE)
lines(seq(1.3, 2.6, by = 0.01), dnorm(seq(1.3, 2.6, by = 0.01), mean = 2, sd = sqrt(tvar)))


#---------------------------
## Interpreting Regression output + Testing
#---------------------------

library(mvtnorm)
set.seed(0)
n<-1000
X<-rmvnorm(n, mean=c(0,0,0), sigma= matrix(c( 1, 0.5, -0.3, 0.5,1, 0, -0.3, 0 , 1), nrow=3) )
beta<-c(0.5,-1.0, 2.0)
x1=X[,1]
x2=X[,2]
x3=X[,3]
y=5+beta[1]*x1+beta[2]*x2 + beta[3]*x3 +rnorm(n)


fitfull<-lm(y ~ x1 + x2 + x3)
s<-summary(fitfull)
s

# Get paramters and associated statistics 
coef<-s$coefficients
coef["x1", "Estimate"]
coef["x1", "Std. Error"]
coef["x1", "t value"]


# Do residual tests
e<-residuals(fitfull)
yhat<-fitted.values(fitfull)

plot(yhat,e)

# Can perform tests with above.

# Single variable tests
#-----------------------
t1<-coef["x1", "t value"] # directly gives the t-value
t2<-coef["x1", "Estimate"]/coef["x1", "Std. Error"] # calculates the t-value

coef["x1", "Pr(>|t|)"] # gives the p-value
(1-pt(t2, df=n-1))*2 # gives the p-value


# Multi variable tests
#-----------------------
s$fstatistic # directly gives the F-statistic value for H0: y= beta_1  HA: Full model

# Can also calculate F-statistics "manually" with the anova function:
fit0<-lm(y ~ 1)
anova(fit0, fitfull)

anova(fitfull,fit0)



# Achtung: For F-tests the two models compared need to be nested!


#---------------------------
# Dummy variables
#---------------------------

library(ISLR)
data(Carseats) #use ?Carseats for an explaination of the dataset
shelveloc=Carseats$ShelveLoc
sales=Carseats$Sales
advertising=Carseats$Advertising

is.factor(shelveloc)
levels(shelveloc)


# fit using automatic coding
fit<-lm(sales~shelveloc+advertising)
summary(fit)

# Interpretation? What is beta1, beta2, beta3?


# Construct dummies -> Watch out for dummy variable trap!
bad<- shelveloc==levels(shelveloc)[1]
medium <- shelveloc==levels(shelveloc)[3]
good<- shelveloc==levels(shelveloc)[2]
a1<-medium*1
a2<-good*1

fit2<-lm(sales~a1+a2+advertising)
summary(fit2)


# How are dummies interpreted?
plot(advertising,sales)
c=coef(fit2)
points(advertising[bad],sales[bad],col="red")
abline(a=c["(Intercept)"],b=c["advertising"],col="red")
points(advertising[medium],sales[medium],col="orange")
abline(a=c["(Intercept)"]+c["a1"],b=c["advertising"],col="orange")
points(advertising[good],sales[good])
abline(a=c["(Intercept)"]+c["a2"],b=c["advertising"])
legend("bottomright", c("bad","medium","good"),col=c("red","orange","black"), lty=1)


# Are dummies needed? Conduct partial F-test
fit0<-lm(sales~advertising)

anova(fit0,fit2)


##################################



#---------------------------
# Transformations
#---------------------------

n<-1000
epsilon<-rnorm(n)
X<-rnorm(n)-5 #rgamma(n,shape=2)
Y<-exp(3+2*X + epsilon) 

plot(X,Y)

plot(X,log(Y))



# Real world example:
airline <- scan("http://stat.ethz.ch/Teaching/Datasets/airline.dat")

# X=t here

plot(airline, type="l")
plot(log(airline), type="l")
t<-1:length(airline)

fit<-lm(log(airline) ~ t)
summary(fit)

et<-residuals(fit)
yhatt=fitted.values(fit)

plot(yhatt,et)
# Important trick: If not sure how the residual plots should look, simply simulate a few for which the assumptions are fullfilled


