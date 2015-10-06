

beta_mixture_fit <- function(x) {
  # note that x is a vector of values whose *distribution* follows a beta mixture distribution
  # do not confuse with x when 1:length(x) vs x shape a beta distribution
  
  library(c(flexmix, betareg))
  d <- data.frame(y=x)
  m <- betamix(y~1|1,data=d,k=1:3)
  mu <- plogis(coef(m)[,1])
  phi <- exp(coef(m)[,2])

  a <- mu*phi
  b <- (1-mu)*phi

  # make a plot of fitting
  # hist(d$y,breaks=0:25/25,freq=F,main="",xlab='y')
  return(list(m=m, mu=mu, phi=phi))
}


kmean <- function(x) {
  m <- kmeans(na.omit(x), centers=c(0.1,0.9))
  # m$cluster
}