#########################################################
# Parametric Bootstrap: 
# Procedure and Model Misspecification
#
# Computational Statistics, Exercise class April 03, 2020
#########################################################

##############################
# Part I: Parametric Bootstrap
##############################

library(boot)

# We consider realizations of a random variable that is 
# Beta(alpha, beta) distributed. 
# We aim to estimate the median as well as the mean.

########### i ###########
# Generate data:
set.seed(8)
n <- 20
alpha <- 2
beta <- 5
xx <- rbeta(n, shape1=alpha, shape2=beta)

# Get an impression of the beta-pdf:
plot(seq(0, 1, length.out = 1000), 
     dbeta(seq(0, 1, length.out = 1000), alpha, beta), 
     xlab = "x",
     ylab = "f(x)", 
     main = sprintf("Beta-density with alpha=%d, beta=%d", alpha, beta))

########### ii ###########
# Maximum likelihood estimates of alpha and beta:
require("MASS")
?fitdistr
(mle.beta <- fitdistr(xx, "beta", start=list(shape1=1, shape2=3))) 
# This outputs parameter estimates (first row)
# and the estimated std-errors (second row). 

# Plot the density with the MLE-estimated parameters on top 
# of the histogram:
par(mfrow=c(1,1))
hist(xx, freq = FALSE, breaks = 50)
lines(x = seq(from = 0, to = max(xx), by = 0.01),
      y = dbeta(x = seq(from = 0, to = max(xx), by = 0.01),
                 shape1 = mle.beta$estimate["shape1"],
                 shape2 = mle.beta$estimate["shape2"]),
      col = 2)

########### iii ###########
# Hard-code parametric bootstrap estimation of the median
# and the mean:
R <- 1000
len.xx <- length(xx)
set.seed(987)
median.star <- rep(NA, R)
mean.star <- rep(NA, R)
for (i in 1:R) {
  x <- rbeta(n = len.xx, 
             shape1 = mle.beta$estimate["shape1"],
             shape2 = mle.beta$estimate["shape2"])
  median.star[i] <- median(x)
  mean.star[i] <- mean(x)
}

# Plot theta*
par(mfrow=c(1,2))
hist(median.star, breaks = 100, freq = FALSE)
abline(v = median(xx), col = 2, lwd = 2)
hist(mean.star, breaks = 100, freq = FALSE)
abline(v = mean(xx), col = 2, lwd = 2)
par(mfrow=c(1,1))

########### iv ###########
# Compute parametric bootstrap estimation of the median
# and the mean using boot:
est.fun <- function(x) {
  return(c(median(x), mean(x)))
}
beta.rg <- function(x, mle) {
  # Function to generate random beta variates.
  # mle will contain the mle of alpha and beta.
  # We simulate from the estimated model.
  rbeta(length(x), shape1 = mle[1], shape2 = mle[2])
}
res.boot <- boot(xx, est.fun, R = R, sim = "parametric", 
                 ran.gen = beta.rg, mle = mle.beta$estimate)
res.boot
res.boot$t0
res.boot$t[1:5,]

########### v ###########
# Confidence intervals: Same as with the ordinary bootstrap!
boot.ci(res.boot, type=c("norm", "basic", "perc"), index=?) #median
boot.ci(res.boot, type=c("norm", "basic", "perc"), index=?) #mean

# To get the bootstrap T confidence intervals, we need to do
# a second level bootstrap:
est.fun <- function(x, index=1:n){
  return(c(median(x[index]), mean(x[index])))
}
est.var.fun <- function(x){           
  mean.med <- est.fun(x)
  boot.sample <- boot(x, est.fun, R=50)
  vars <- c(var(boot.sample$t[,1]), var(boot.sample$t[,2]))
  return(c(mean.med, vars))
}
beta.rg <- function(x, mle) {
  # Function to generate random beta variates.
  # mle will contain the mle of alpha and beta.
  rbeta(length(x), shape1 = mle[1], shape2 = mle[2])
}
res.boot2 <- boot(xx, est.var.fun, R = R, sim = "parametric", 
                  ran.gen = beta.rg, mle = mle.beta$estimate)

res.boot2
res.boot2$t0
res.boot2$t[1:5,]

boot.ci(res.boot2, type="stud", index=?) #median
boot.ci(res.boot2, type="stud", index=?) #mean


####################################
# Part II: Model Misspecification in 
#          the Parametric Bootstrap
###################################

########### i ###########
# Generate data:
set.seed(8)
n <- 30
alpha <- 1
beta <- 4
xx <- rbeta(n, shape1=alpha, shape2=beta)

# Get an impression of the beta-pdf:
plot(seq(0, 1, length.out = 1000), 
     dbeta(seq(0, 1, length.out = 1000), alpha, beta), 
     xlab = "x",
     ylab = "f(x)", 
     main = sprintf("Beta-pdf with alpha=%d, beta=%d", alpha, beta))

########### ii ###########
# Maximum likelihood estimates of alpha and beta:
require("MASS")
?fitdistr
(mle.lnorm <- fitdistr(xx, "log-normal")) 
# This outputs parameter estimates (first row)
# and the estimated std-errors (second row)

# Get an impression of the log-normal-pdf:
plot(seq(0, 1, length.out = 1000), 
     dlnorm(seq(0, 1, length.out = 1000), 
            meanlog = mle.lnorm$estimate["meanlog"], 
            sdlog = mle.lnorm$estimate["sdlog"]), 
     xlab = "x",
     ylab = "f(x)", 
     main = sprintf("Lognormal-pdf with mu=%f, sigma=%f", 
                    mle.lnorm$estimate["meanlog"], 
                    mle.lnorm$estimate["sdlog"]))

# Plot the density with the MLE-estimated parameters on top 
# of the histogram:                                          
hist(xx, freq = FALSE, breaks = 50)
lines(x = seq(from = 0, to = max(xx), by = 0.01),
      y = dbeta(x = seq(from = 0, to = max(xx), by = 0.01),
                shape1 = alpha,
                shape2 = beta),
      col = 2)
lines(x = seq(from = 0, to = max(xx), by = 0.01),
      y = dlnorm(x = seq(from = 0, to = max(xx), by = 0.01),
                 meanlog = mle.lnorm$estimate["meanlog"], 
                 sdlog = mle.lnorm$estimate["sdlog"]),
      col = 4)
legend("topright", legend=c("beta", "lognormal"), col=c(2,4), lty=1)

########### iii ###########
# Parametric bootstrap:
R <- 10000
len.xx <- length(xx)
set.seed(987)
param.quant.star <- rep(NA, R)
param.mean.star <- rep(NA, R)
for (i in 1:R) {
  x <- rlnorm(n = len.xx, 
             meanlog = mle.lnorm$estimate["meanlog"],
             sdlog = mle.lnorm$estimate["sdlog"])
  param.quant.star[i] <- quantile(x, probs = 0.95)
  param.mean.star[i] <- mean(x)
}

########### iv ###########
# non-parametric bootstrap
nonparam.quant.star = rep(NA, R)
nonparam.mean.star = rep(NA, R)
for(i in 1:R){
  bootstrap_sample <- xx[sample(1:n, n, replace=TRUE)]

  nonparam.quant.star[i] = quantile(bootstrap_sample, probs=0.95)
  nonparam.mean.star[i] = mean(bootstrap_sample)
}

########### v ###########
# Investigate estimated standard errors:
sim.beta.quant.mean <- function(){
  x <- rbeta(n, shape1=alpha, shape2=beta)
  return(c(quantile(x, probs=0.95), mean(x)))
}
reps <- 10000
sim.result <- replicate(reps,sim.beta.quant.mean()) #(2 x reps)-matrix where first
                                                    #row contains the reps realizations
                                                    #of the quant and the second row
                                                    #contains the means
rownames(sim.result) <- c("95%-quant", "mean")
sim.result[,1:5]

true.sd.quant <- sd(sim.result[1,])
true.sd.mean <- sd(sim.result[2,])
param.sd.quant <- sd(param.quant.star)
param.sd.mean <- sd(param.mean.star)
nonparam.sd.quant <- sd(nonparam.quant.star)
nonparam.sd.mean <- sd(nonparam.mean.star)

wrapup <- rbind(c(true.sd.quant, param.sd.quant, nonparam.sd.quant), 
                 c(true.sd.mean, param.sd.mean, nonparam.sd.mean))
rownames(wrapup) <- c("0.95 quantile", "mean")
colnames(wrapup) <- c("true", "parametric", "nonparametric")
wrapup
