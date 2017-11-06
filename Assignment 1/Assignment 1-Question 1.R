######################################
###COMM 590  Assignment 1 Question 1
###Qiyuan Wang
######################################
library(ggplot2)
library(expm)
library(matlib)

###multivariate normal distribution
DistNorm=function(n,mean,Var_Cov){
  #this function generate n*m matrix with each row follows m dimensional multivariate normal distribution
  #n is the number of observation
  #m is random variable dimension
  #mean is the mean vector
  #Var_Cov is variance-covariance matrix for the ramdom variable
  m=length(mean)
  x=rnorm(n*m,0,1)
  xmulti=matrix(x,ncol=m)
  xmulti=xmulti%*%sqrtm(Var_Cov)+t(replicate(n,mean))
  return(xmulti)
}
#show the results
Var_Cov=matrix(c(1,-1,-1,2),nrow=2)
x=DistNorm(10000,c(1,2),Var_Cov)
qplot(x[,1],x[,2], geom='bin2d')+
  labs(x='x',y='y')

###gamma and inverse gamma distribution
#gamma distrubtion
#reference:https://www.statlect.com/probability-distributions/gamma-distribution
DistGamma=function(k,theta){
  #this function generates random variable from gamma distribution
  #It should be noted that this function only work for gamma distribution with theta being integer
  x=rnorm(theta,0,1);
  x=sum(x^2);
  xgamma=(k/theta)*x
  return(xgamma)
}
#show the results
test=numeric(10000)
for(i in 1:10000){
  test[i]=DistGamma(0.5,7)
}
qplot(test,fill=I('steelblue'))+
  labs(x='Gamma Distribution')

#inverse gamma
DistInverseGamma=function(k,theta){
  #this function generates random variable from inverse gamma distribution
  xinversegamma=DistGamma(k,theta)
  return(1/xinversegamma)
}
#show the results
test=numeric(10000)
for(i in 1:10000){
  test[i]=DistInverseGamma(0.5,7)
}
qplot(test,fill=I('steelblue'),binwidth=0.3)+
  labs(x='Inverse Gamma Distribution')

###Whishart and inverse whishart distribution
DistWhishart=function(n,Var_Cov){
  #this function generates random matrix that follows whishart distribution
  #it uses the multivariate random normal variable generating function DistNorm
  m=nrow(Var_Cov)
  #generate the random variable
  xmulti=DistNorm(n,rep(0,times=m),Var_Cov)
  xwhishart=t(xmulti)%*%xmulti
  return(xwhishart)
}
#show the results
DistWhishart(3,Var_Cov)

#inverse Whishart distribution
DistInverseWhishart=function(n,Var_Cov){
  #this function generates random matrix that follows inverse whishart distribution
  #it uses whishart generating function DistNorm
  xinversewhishart=inv(DistWhishart(n,Var_Cov))
  return(xinversewhishart)
}
#show the results
DistInverseWhishart(3,Var_Cov)

