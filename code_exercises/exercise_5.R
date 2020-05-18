n = 10000

probability.bootstrap <-function(n){
  y <- 1 - (1 - 1/n)^n
}
data <- 1:n
plot(data, probability.bootstrap(data), type = 'l', log = "x")

