######################################
###COMM 590  Assignment 1 Question 2
###Qiyuan Wang
######################################

###I made some new changes
library(invgamma)
library(MASS)
library(ggplot2)
library(Matrix)
library(MCMCpack)
library(reshape2)

###############################################
#Univariate Regression Model
###############################################
univariate=function(X,Y,N,mean,variance,alhpa,beta){
  #This function use Gibbs Sampler to draw from posterior distribution for linear models
  #X,Y are idependent variables matrix and Y is dependent variable vector
  #N is the number of iteration
  #mean and variance are the hyper parameters for regression coefficients' normal prior
  #alpha and beta are the hyper parameters for standard error' gamma prior
  
  nparams=ncol(X)+2#number of parameters to be estimated=number of DV+intercept+variance
  nobs=nrow(X)#number of observations
  X=cbind(matrix(1,nrow=nobs,ncol=1),X)
  
  #set the matrix to store draws
  draw=matrix(0,nrow=N,ncol=nparams)
  
  #set the starting value for draws
  draw[1,]=1
  

  betahat=solve(t(X)%*%X)%*%t(X)%*%Y
  #Gibbs Sampler
  for(i in 2:N){
    #Firstly, draw from gamma
    #the estimate for OLS
    draw[i,1]=rinvgamma(1, shape = (nobs/2 + alhpa), 
                        scale =.5*sum((Y-(X%*%betahat))^2)+beta)
    #Secondly, draw from normal distribution
    #calculate the mean and variance for posterior distribution
    meannew=solve(t(X)%*%X+solve(variance/draw[i,1]))%*%(t(X)%*%Y+solve(variance/draw[i,1])%*%mean)
    variancenew=draw[i,1]*solve(t(X)%*%X+solve(variance/draw[i,1]))
    #draw from multivariate normal distribution
    draw[i,2:nparams]=mvrnorm(1,meannew,variancenew)
    
    #show the iteration results
    cat(paste("Iteration",i,sep=' '),":",draw[i,], "\n")
  }
  return(draw)
}

########test the function
#generate some data
X=as.matrix(rnorm(1000,0,1))
epsilon=as.matrix(rnorm(1000,0,2))
Y=2+2*X+epsilon
#set up hyperparamters
mean=c(1,1)
variance=matrix(c(1,0,0,1),nrow=2)
alhpa=1
beta=1
#iteration times
N=1000
#run the function
test=univariate(X,Y,N,mean,variance,alhpa,beta)
test=as.data.frame(test)
test$Iteration=1:1000
names(test)=c('Sigma','Beta0','Beta1','Iteration')
test=melt(test,id.vars = 4)
names(test)=c('Iteration','Variable','Value')
ggplot(data=test[test$Iteration>=200,],aes(x=Iteration,y=Value,color=Variable))+
  geom_line()

###############################################
#Multivariate Regression Model
###############################################
multivariate=function(X,Y,N,nu,v,beta,sigma,A){
  #This function use Gibbs Sampler to draw from posterior distribution for linear models
  #X,Y are idependent variables matrix and Y is dependent variable vector
  #N is the number of iteration
  #mean and variance are the hyper parameters for regression coefficients' normal prior
  #alpha and beta are the hyper parameters for standard error' gamma prior
  
  ny=ncol(Y)#number dependent variables
  nx=ncol(X)+1#number independent variables
  nobs=nrow(Y)#number of observations
  
  X=cbind(matrix(1,nrow=nobs,ncol=1),X)
  
  #set the matrix to store draws
  drawsigma=matrix(0,nrow=N*ny,ncol=ny)#the draw for variance
  drawbeta=matrix(0,nrow=N*ny*nx)#the draw for beta
  
  #set the starting value for draws
  drawsigma[1:ny,1:ny]=1
  drawbeta[1:(ny*nx),]=1
  
  #calculate relevant matrix
  Btilda=solve(t(X)%*%X+A)%*%(t(X)%*%Y+A%*%beta)
  S=t(Y-X%*%Btilda)%*%(Y-X%*%Btilda)+t(Btilda-beta)%*%A%*%(Btilda-beta)
  #Gibbs Sampler
  for(i in 2:N){
    #Firstly, draw from invert whishart distribution
    tempsigma=riwish(nu+nobs,v+S)
    drawsigma[((i-1)*ny+1):(i*ny),]=tempsigma
    #Secondly, draw from normal distribution
    #calculate the mean and variance for posterior distribution
    meannew=c(Btilda)
    variancenew=kronecker(tempsigma, solve(t(X)%*%X+A), FUN = "*")
    #draw from multivariate normal distribution
    drawbeta[((i-1)*(ny*nx)+1):(i*(ny*nx)),]=mvrnorm(1,meannew,variancenew)
    
    #show the iteration results
    cat("Iteration",i,"for Sigma:",tempsigma, "\n")
    cat("Iteration",i,"for Beta:",drawbeta[((i-1)*(ny*nx)+1):(i*(ny*nx)),], "\n")
  }
  return(list(drawsigma,drawbeta))
}

########test the function
#generate some data
X=matrix(rnorm(2000,0,1),nrow=1000)
beta=matrix(c(1,2,3,4,5,6),nrow=3)
epsilon=matrix(rnorm(2000,0,1),nrow=1000)
Y=cbind(matrix(1,nrow=1000,ncol=1),X)%*%beta+epsilon
#set up hyperparamters
nu=2;
v=diag(rep(1,2))
beta=matrix(1,nrow=3,ncol=2)
sigma=diag(rep(1,2))
A=matrix(1,nrow=3,ncol=3)
#iteration times
N=1000
#run the function
test=multivariate(X,Y,N,nu,v,beta,sigma,A)
test1=as.data.frame(test[1])
test2=as.data.frame(test[2])
temp=data.frame(Iteration=1:1000,
                sigma1=test1[seq(1,1999,by=2),1],
                Sigma2=test1[seq(1,1999,by=2),2],
                Sigma3=test1[seq(2,2000,by=2),1],
                Sigma4=test1[seq(2,2000,by=2),2],
                Beta1=test2[seq(1,5999,by=6),],
                Beta2=test2[seq(2,5999,by=6),],
                Beta3=test2[seq(3,5999,by=6),],
                Beta4=test2[seq(4,5999,by=6),],
                Beta5=test2[seq(5,5999,by=6),],
                Beta6=test2[seq(6,6000,by=6),])
temp=melt(temp,id.vars = 1)
names(temp)=c('Iteration','Variable','Value')
ggplot(data=temp[temp$Iteration>=200,],aes(x=Iteration,y=Value,color=Variable))+
  geom_line()
