#########################################################
# Estimating the Accuracy of a Linear Regression Model
# Simulate from model with heteroscedastic errors 
# (non-constant variance)
#
# Computational Statistics, Exercise class March 27, 2020
#########################################################

library(boot)


#########################################
# Part I: Problem statement
#########################################
set.seed(999)
n <- 100 #sample size

# Function that returns a generated data set (x,y):
generate.data <- function(){
  x <- runif(n, 0, 100)
  ?runif
  
  eps <- rnorm(n,0,abs(x-50))
  y <- 50 + 0.3*x + eps
  return(data.frame(cbind(x,y)))
}
# Our data set:
dat <- generate.data()

# Fit a linear regression:
fit1 <- lm(y~x, data=dat)
summary(fit1) 

 # Visualize our data:
par(mfrow=c(1,1))
plot(y~x,data=dat)
abline(50, 0.3) #truth
abline(coef(fit1), col='red') #estimate


#########################################
# Part II: Bootstrap estimates
#########################################

# Function used for bootstrap needs to be of this format.
# Tt accepts the array of indices which is used to 
# sample the data:
boot.fn <- function(data,index){
  return(coef(lm(y~x, data=data, subset=index)))
}

# Compute Bootstrap estimates of the coefficients:
?boot
boot.fit <- boot(dat, statistic=boot.fn, R=5000)
boot.fit$t0 #original estimate (using all data)
boot.fit$t[1:5,] #the outputs of statistic function
boot.fit

# Explanation of the output of boot.fit:
# -) "$t0": This contains the values of our statistic(s) in 
#    the original, full data set.
# -) "$t": This contains the R (=B) values of our statistic(s)
#    that were generated in the individual bootstrap repetitions.
#    That is, $t contains the \theta^(*i)'s, i=1,...,R.
# -) "t1* & t2*": In our case, the t1* abbreviates the intercept
#    and the t2* the slope.
# -) "original": This is the same as $t0.
# -) "bias": For instance, the value given in t1* equals
#    mean(boot.fit$t[,1])-boot.fit$t0[1]
# -) "std. error": This is the standard error/standard deviation
#    of the bootstrap realizations. For instance, the value given 
#    in t1* equals sd(boot.fit$t[,1]).


#########################################
# Part III: Bootstrap confidence 
#           intervals
#########################################

confint(fit1) #confidence interval from lm
?boot.ci
 boot.ci(boot.fit, conf=0.90, type="basic") #bootstrap CI, conf level is 
                                            #by default 0.95
boot.ci(boot.fit, index=c(2), type="basic") #index needed for other variables 
                                            #in output of boot

(ci.basic <- boot.ci(boot.fit, index=c(1), type="basic"))
(ci.perc <- boot.ci(boot.fit, index=c(1), type="perc"))
(ci.norm <- boot.ci(boot.fit, index=c(1), type="norm"))
boot.ci(boot.fit, index=c(1), type="stud") #error, need to give variances 
                                           #for t statistics

# For "stud", our bootstrap function needs to return 
# an estimate of the variances as well:
coef.fn <- function(data, index){
  return(coef(lm(y~x,data=data,subset=index)))
}
boot.fn <- function(data,index){
  coefs <- coef.fn(data, index)
  boot.sample <- boot(data[index, ], coef.fn, R=50)
  vars <- c(var(boot.sample$t[,1]), var(boot.sample$t[,2]))
  return(c(coefs, vars))
}
boot.fn(dat, 1:100)

boot.fit2 <- boot(dat, boot.fn, R=300) #takes time because we have a bootstrap 
                                       #inside a bootstrap
boot.fit2$t0 
boot.fit2$t[1:5,]

# index[1] tells us in which position of boot.fit2$t0 
# is the estimator is and index[2] tells us where is its variance is. 
# if our bootstrap statistic returns two values and we do not specify index, 
# then by default first comes the estimator and then its variance.
(ci.stud <- boot.ci(boot.fit2, index=c(1, 3), type="stud")) #for intercept
boot.ci(boot.fit2, index=c(2, 4), type="stud") #for slope
# Equivalent way of computing the confidence intervals:
boot.ci(boot.fit2, type="stud", 
        var.t0=boot.fit2$t0[3], 
        var.t=boot.fit2$t[,3]) #for intercept
boot.ci(boot.fit2, type="stud", 
        var.t0=boot.fit2$t0[4], 
        var.t=boot.fit2$t[,4]) #for slope

# Compute all CIs at the same time:
(ci.all <- boot.ci(boot.fit2, index=c(1, 3), 
                   type=c("basic", "perc", "norm", "stud")))
ci.all[["normal"]][2:3]
ci.all[["percent"]][4:5]
ci.all[["basic"]][4:5]
ci.all[["student"]][4:5]


#########################################
# Part IV: Coverage of the Bootstrap
#          confidence intervals for the 
#          intercept
#########################################

true.par <- 50 #the true intercept
ci.basic
ci.perc
ci.norm
ci.stud

# Repeate the above procedure of computing the intervals
# ci.basic, ci.perc, ci.norm and ci.stud a total number
# of n.sim times. Then look at proportion of the intervals
# where true.param is to the left and where true.param is
# to the right.
# We then compute the 4 CIs a total of n.sim times for
# different values of n. 
# Please see the slides for the plots.


#########################################
# Part V: Compute the CIs by hand
#########################################

R <- length(boot.fit2$t[,1])
theta.hat <- boot.fit2$t0[1]
theta.star.bar <- mean(boot.fit2$t[,1])

# Bootstrap T CI: 
sd.hat <- sqrt(sum((boot.fit2$t[,1]-theta.star.bar)^2)/(R-1))
t.stats <- (boot.fit2$t[,1]-theta.hat)/sqrt(boot.fit2$t[,3])
quant <- quantile(t.stats, probs = c(0.975, 0.025))
theta.hat - sd.hat*quant

boot.ci(boot.fit2, index=c(1, 3), type = "stud")[["student"]][4:5]


# Normal CI with bias correction:
z.quant <- qnorm(0.975, 0, 1)
sd.hat <- sqrt(sum((boot.fit2$t[,1]-theta.star.bar)^2)/(R-1))
c(2*theta.hat - theta.star.bar - z.quant*sd.hat, 
  2*theta.hat - theta.star.bar + z.quant*sd.hat)

boot.ci(boot.fit2, index=c(1), type="norm")[["normal"]][2:3]


# Percentile (=Quantile) CI:
quantile(boot.fit2$t[,1], probs = c(.025, .975))
boot.ci(boot.fit2, index=c(1), type="perc")[["percent"]][4:5]


# Basic (=Reversed Quantile)
quant <- quantile(boot.fit2$t[,1]-theta.hat, 
                  probs = c(0.975, 0.025))
theta.hat-quant

boot.ci(boot.fit2, index=c(1), type = "basic")[["basic"]][4:5]
