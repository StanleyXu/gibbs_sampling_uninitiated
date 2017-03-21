#### GIBBS SAMPLING FOR THE UNINITIATED

source("gibbs.R")

################################################################################
#### Set up
################################################################################
# number of documents in the corpus
N=100
# number of words in the vocabulary
V=5

################################################################################
#### Simulation: Generate truth
################################################################################
# hyperpara of the beta distribution
#gamma.pi.1=2
#gamma.pi.0=8
# hyperpara vector for the multinomial prior
#gamma.theta=rep(1,V)

truth=simu(N=N, V=V, gamma.pi.1=6, gamma.pi.0=4, gamma.theta=rep(1,V))

################################################################################
#### Gibbs Sampling
################################################################################

result=gibbs(N, V, W=truth$W, iter=20)

################################################################################
#### Compare the sampling result with truth
################################################################################

## label flip
if(sum(abs(result$theta[1,] - truth$theta[1,])) > sum(abs(result$theta[1,] - truth$theta[2,]))){
  result$theta=result$theta[2:1,]
  tt=result$per.iter$theta.0
  result$per.iter$theta.0=result$per.iter$theta.1
  result$per.iter$theta.1=tt; rm(tt)
}

print(result$theta)
print(truth$theta)

mse.0=apply(result$per.iter$theta.0, 1, FUN=function(x) mean((x-truth$theta[1,])^2))
mse.1=apply(result$per.iter$theta.1, 1, FUN=function(x) mean((x-truth$theta[2,])^2))

plot(mse.0, ylim=range(c(mse.0, mse.1)), type='l', col='blue', ylab="MSE", xlab="iteration")
points(mse.1, type='l', col='red')
legend("topright", legend=c("theta_0", "theta_1"), col=c("blue", "red"), lty=1, bty='n')
