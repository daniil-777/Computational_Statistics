set.seed(0)
n<-100
z1<-rnorm(n)
z2<-rnorm(n)
M=matrix(c(1,1,0.1,-0.1),2,2)
X=t(M%*%rbind(z1,z2))
beta<-c(0.5,-1.0)
x1=X[,1]
x2=X[,2]
y=5+beta[1]*x1+beta[2]*x2 +rnorm(n)

n <- 100 #the amoubt of data in x1/x2 (X)
p <- 3 #intercept and two variables

#claculation of esimation parameters manually - 1-st way
fit1<-lm(y~x1+x2)
######################################################
#first way to calculate estimated_std
summary(fit1)
s<-summary(fit1)
coef<-s$coefficients
intercept <- coef["(Intercept)", "Estimate"]
coef_1 <- coef["x1", "Estimate"]
coef_2 <- coef["x2", "Estimate"]
coef_vector <- t(cbind(intercept, coef_1, coef_2)) # betta hat
y_pred = X %*% coef_vector #y_hat
estimated_std <- 1/(n - p)*sum((y - y_pred)^2) #std.hat^2
########################################################
#second way to calculate estimated_std

X <- cbind(1,x1, x2) #design matrix
X <-as.matrix(cbind(1,x1, x2))
XtX.inv <- solve(t(X)%*%X)
betta.hat <- XtX.inv %*% t(X) %*% y

y.hat <- X %*% betta.hat
res = y - y.hat
estimated_std.betta1 <- 1/(n - p)*sum(res^2) #std "overall"
se.betta_1 = sqrt(estimated_std.betta1*XtX.inv[2,2])
t_value.betta_1 <- betta.hat[2]/se.betta_1

#standart error calculation by just calling and manually - se(betta_{k}.hat)
tvar <- coef["x1", "Std. Error"] 
tvar_manual <-sqrt(estimated_std*XtX.inv[2,2]) 


#p value calculation
?pt
2*pt(abs(t_value.betta_1), df=n-p, lower=FALSE)


(1-pt(t_value.betta_1, df=n-1))*2 # gives the p-value


