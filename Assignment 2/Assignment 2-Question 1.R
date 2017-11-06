######################################
###COMM 590  Assignment 2 Question 1
###Qiyuan Wang
######################################

# Library for sampling from Multivariate Normal distribution
require(mvtnorm)
# Library for sampling from Truncated Normal distribution
require(truncnorm)


###Binary Probit Model
BiProbit=function(X,Y,N,mean,variance){
  #This function use Gibbs Sampler to draw from posterior distribution for linear models
  #X,Y are idependent variables matrix and Y is dependent variable vector
  #N is the number of iteration
  #mean and variance are the hyper parameters for regression coefficients' normal prior
  #alpha and beta are the hyper parameters for standard error' gamma prior
  
  nparams=ncol(X)+1#number of parameters to be estimated=number of DV+intercept+latent variable
  nobs=nrow(X)#number of observations
  X=cbind(matrix(1,nrow=nobs,ncol=1),X)
  
  #set the matrix to store draws for beta
  draw=matrix(0,nrow=N,ncol=nparams)
  #set the matrix tos tore draws for latent variable
  latent=matrix(0,nrow=nobs,ncol=1)
  
  #set the starting value for draws
  draw[1,]=1
  
  
  #Gibbs Sampler
  for(i in 2:N){
    #Firstly, draw for latent variable from truncated normal distribution
    #the mean for truncated normal distribution
    meantru=X%*%draw[i-1,1:nparams]
    #the estimate for OLS
    latent[which(Y==1),1]=rtruncnorm(1, a=0, b=Inf, mean = meantru[which(Y==1),1], sd = 1)
    latent[which(Y==0),1]=rtruncnorm(1, a=-Inf, b=0, mean = meantru[which(Y==0),1], sd = 1)
    
    #Secondly, draw for beta from normal distribution
    #calculate the mean and variance for posterior distribution
    meantheta=solve(t(X)%*%X+solve(variance))%*%(t(X)%*%latent+solve(variance)%*%mean)
    variancetheta=solve(t(X)%*%X+solve(variance))
    #draw from multivariate normal distribution
    draw[i,]=mvrnorm(1,meantheta,variancetheta)
    
    #show the iteration results
    cat(paste("Iteration",i,sep=' '),":",draw[i,], "\n")
  }
  return(draw)
}


##############test the model
#generate data
beta=c(-1,0.5,1.5)
#generate X
X=rnorm(2000)
X=matrix(X,nrow=1000)
Xtemp=cbind(rep(1,times=1000),X)
#chocie probability
PChoice=pnorm(Xtemp%*%beta)
temp=runif(1000)
Y=as.matrix(as.numeric(temp<=PChoice))

#run the model
N=10000
mean=c(1,1,1)
variance=diag(c(1,1,1))
test=BiProbit(X,Y,N,mean,variance)
qplot(1:10000,test[,1])
colMeans(test[1000:10000,])
