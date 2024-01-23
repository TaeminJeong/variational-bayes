set.seed(1)
N <- 50
y <- rnorm(n=N, mean=10, sd=sqrt(5))
mean_y <- mean(y); mean_y

#### VB(MFVB) ####
## setting initial value ##
mu_mu0 <- mu_mu <- 0
sigma2_mu0 <- sigma2_mu <- 1
alpha0 <- alpha <- 1
beta0 <- beta <- 1

## running VB ##
for(i in 1:10000){
  alpha <- alpha0 + N/2
  beta <- beta0 + 1/2*sum(y^2) - N*mean_y*mu_mu + N/2*(mu_mu^2+sigma2_mu)
  mu_mu <- (mu_mu0/sigma2_mu0 + N*mean_y*alpha/beta)/(1/sigma2_mu0+N*alpha/beta)
  sigma2_mu <- (1/sigma2_mu0+N*alpha/beta)^(-1)
}

#### MCMC ####
## model for MCMC ##
library(nimble)
code <- nimbleCode({
  for(i in 1:N){
    y[i] ~ dnorm(mean=mu, var=sigma2)
  }
  mu ~ dnorm(mean=0, sd=1)
  sigma2 ~ dinvgamma(shape=1,scale=1)
})
constants = list(N=N)
data = list(y=y)
monitors = c('mu','sigma2')
inits <- list(mu=0, sigma2=1)

## running MCMC ##
params <- nimbleMCMC(code=code, constants=constants,inits=inits,data=data,
                     monitors=monitors,niter=10000,nburnin=1000)

## result ##
plot(density(params[,1]),xlim=c(8,11),ylim=c(0,1.5),main='density of mu',ylab='',xlab='')
lines(seq(8,11,length=100),
      dnorm(seq(8,11,length=100),
            mean = mu_mu,sd = sqrt(sigma2_mu)),
      xlim=c(8,11),ylim=c(0,1.5),col='red')
legend('topright',legend=c('MCMC','VB'),lty=1,col=c('black','red'))

plot(density(params[,2]),xlim=c(2,8),ylim=c(0,0.55),main='density of sigma2',ylab='',xlab='')
lines(seq(2,8,length=100),
      dinvgamma(seq(2,8,length=100),
                shape = alpha,scale = beta),
      xlim=c(2,8),ylim=c(0,0.55),col='red')
legend('topright',legend=c('MCMC','VB'),lty=1,col=c('black','red'))
