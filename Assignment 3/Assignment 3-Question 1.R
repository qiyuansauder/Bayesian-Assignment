#load library
library(mvtnorm)
library(data.table)
library(invgamma)
library(MCMCpack)

####read in data
process=read.csv('D:/Dropbox/1 Courses/11-So Eun Chunhua-Structural Model/Chunhua Assignment/googleCLV_data/process.csv')
amount=read.csv('D:/Dropbox/1 Courses/11-So Eun Chunhua-Structural Model/Chunhua Assignment/googleCLV_data/amount.csv')
customer=read.csv('D:/Dropbox/1 Courses/11-So Eun Chunhua-Structural Model/Chunhua Assignment/googleCLV_data/customer.csv')

###the individual likelihood for mu and lamda
likelihood1=function(mu,lamda,x,tx,t){
  ll=((lamda^x*mu)/(lamda+mu))*exp(-(lamda+mu)*tx)+
    (lamda^(x+1))/(mu+lamda)*exp(-(lamda+mu)*t)
  return(ll)
}
###the coditional density for mu and lamda given b,G and sigma
likelihood2=function(mu,lamda,mean,variance){
  temp=cbind(mu,lamda,mean)
  p=apply(temp,1,function(x)dmvnorm(x[1:2],mean=x[3:4],sigma=variance))
  return(p)
}

###the likelihood ratio for transaction rate
likelihood3=function(mu1,lamda1,mu2,lamda2,x,tx,t){
  numerator=(((lamda1^x)*mu1)/(lamda1+mu1))*exp(-(lamda1+mu1-mu2-lamda2)*tx)+
    (lamda1^(x+1))/(mu1+lamda1)*exp(-(mu1+lamda1)*t+(mu2+lamda2)*tx)
  denominator=(((lamda2^x)*mu2)/(lamda2+mu2))+
    (lamda2^(x+1))/(mu2+lamda2)*exp(-(mu2+lamda2)*(t-tx))
  return(numerator/denominator)
}

###multiply the two likelihood
likelihood=function(mu,lamda,mean,variance,x,tx,t){
  p=likelihood1(mu,lamda,x,tx,t)*likelihood2(mu,lamda,mean,variance)
  return(p)
}

###the iteration process
draw=function(iteration,beta0,B0,shape,scale,nu,v,G0,A,amount,customer,process){
  ###input
  #iteration is the number of draws
  #beta0,B0 is the mean and variance of normal prior for beta
  #shape,scale is the shape and scale parameter of gamma prior for epsilon's variance
  #nu,v is the parameter of inverse whishart prior for sigma
  #G0,A is the parameter of normal prior for G
  #amount,customer, process is the datasets
  
  ###extract data from original dataset
    #customer dataset:customer attributes
    X=customer[,2:5]
    X=as.matrix(cbind(matrix(1,nrow=408,ncol=1),X))
    #amount datasset: customer expenditure
    z=amount[,c(1,3)]
    z$value=log(z$amount)
    z=as.data.table(z)
    d=amount$t
    ind=which(amount$t>1)
    #process dataset: customer purchase history
    x=process$n#total transactions
    tx=process$tx#last time purchase
    t=process$T#total length of observation
    #number of customers
    n=length(process$custID)
  
  ###create matrix to store the draw
    #parameters for hierachicial parameters
    tempparam=matrix(0,nrow=iteration,ncol=17)
    #variance covariance matrix
    tempvariance=matrix(0,nrow=iteration*3,ncol=3)
  
  ###set the starting value
    #paramters to store
    tempparam[1,]=1
    tempvariance[1:3,1:3]=diag(3)
    #latent variables for each individual,including log mu,log lamda and b
    eta=matrix(1,nrow=n,ncol=2)
    b=matrix(1,nrow=n,ncol=1)
  

  ###drawing values from posterior distribution
    for(i in 2:iteration){
      #extract prior draw at i-1
      G=as.matrix(rbind(tempparam[i-1,1:5],tempparam[i-1,6:10],tempparam[i-1,11:15]))
      Sigma=tempvariance[((i-2)*3+1):((i-1)*3),]
      beta=tempparam[i-1,16]
      sigmaepsilon=tempparam[i-1,17]
      
      #draw for mu and lamda
        #decompose the variance covariance matrix
        Sigma11=Sigma[1:2,1:2]
        Sigma12=as.matrix(Sigma[1:2,3])
        Sigma21=t(as.matrix(Sigma[3,1:2]))
        Sigmab=Sigma[3,3]
        #calculate the mean value for mu,lamda and b
        thetabar=X%*%t(G)
        #get the conditional mean and variance for mu and lamda
        etatilda=eta+(b-thetabar[,3])%*%t(Sigma12)*(1/Sigmab)
        Sigma11tilda=Sigma11-Sigma12%*%Sigma21*(1/Sigmab)
        #using metropolis-hasting to draw mu and lamda
          #generate steps
          etaprime=rnorm(n*2,mean=0,sd=0.05)
          etatemp=eta+cbind(etaprime[1:n],etaprime[(n+1):(n*2)])
          #calculate the acceptance probability
          p1=likelihood3(exp(etatemp[,1]),exp(etatemp[,2]),exp(eta[,1]),exp(eta[,2]),x,tx,t)
          p2=likelihood2(exp(etatemp[,1]),exp(etatemp[,2]),etatilda,Sigma11tilda)/
            likelihood2(exp(eta[,1]),exp(eta[,2]),etatilda,Sigma11tilda)
          p=p1*p2
          p=cbind(matrix(1,nrow=408,ncol=1),p)
          paccept=apply(p,1,min)
          #use uniformly distributed random number to decide whether accept the proposal
          accept=paccept>=runif(408)
          eta[which(accept==T),]=etatemp[which(accept==T),]
          
          
      #draw for b
        #calculate the mean and variance for posterior
        sigmabhat=Sigmab-Sigma21%*%solve(Sigma11)%*%Sigma12
        bhat=thetabar[,3]+(eta-thetabar[,1:2])%*%t(Sigma21)*(1/Sigmab)
        sigmabtilda=1/(1/sigmabhat+x*(1/sigmaepsilon))
        #sum over z value
        z$value[ind]=log(z$amount[ind])-beta*log(d[ind])
        zvalue=z[,.(value=sum(value)),by=.(custID)]
        btilda=sigmabtilda*((1/sigmabhat[1,1])*bhat+(1/sigmaepsilon)*(zvalue$value))
        btemp=cbind(btilda,sigmabtilda)
        b=apply(btemp,1,function(x)rnorm(1,mean=x[1],sd=x[2]))
     
      #draw for beta and delta
      zy=log(z$amount)-b[amount$custID]
      zx=log(d)
      betahat=t(zx)%*%zy/(t(zx)%*%zx)
      sigmaepsilon=rinvgamma(1, shape = (length(zx)/2 + shape), 
                             scale =.5*t((zy-(zx%*%betahat)))%*%(zy-(zx%*%betahat))/sigmaepsilon+scale)
      meannew=(1/(t(zx)%*%zx/sigmaepsilon+1/B0))%*%(t(zx)%*%zy/sigmaepsilon+(1/B0)%*%beta0)
      variancenew=sigmaepsilon*(1/(t(zx)%*%zx+1/(B0/sigmaepsilon)))
      #draw from multivariate normal distribution
      beta=rnorm(1,meannew,variancenew)
  
      #draw for G and sigma 
      Y=cbind(eta,b)
      RA=chol(A)
      W=rbind(X,RA)
      Z=rbind(Y,RA%*%G0)
      Btilda=solve(t(X)%*%X+A)%*%(t(X)%*%Y+A%*%G0)
      S=t(Z-W%*%Btilda)%*%(Z-W%*%Btilda)
      Sigma=rWishart(1,nu+n,solve(v+S))
      Sigma=solve(Sigma[1:3,1:3,1])
      meannew=c(Btilda)
      variancenew=kronecker(Sigma, solve(t(X)%*%X+A), FUN = "*")
      #draw from multivariate normal distribution
      G=mvrnorm(1,meannew,variancenew)
      
      ###store the results
      tempparam[i,1:15]=G
      tempparam[i,16]=beta
      tempparam[i,17]=sigmaepsilon
      tempvariance[((i-1)*3+1):(i*3),1:3]=Sigma

      ###print the results
        #every 100 iteration print the results
        if(i%%100==0){
          cat("Iteration",i,": ",tempparam[i,],c(Sigma), "\n")
        }
    }
    return(list(tempparam,tempvariance))
}


###test the code
iteration=10000
shape=1
scale=1
beta0=1
B0=1
A=diag(5)
G0=matrix(2,nrow=5,ncol=3)
nu=2
v=diag(3)
test1=draw(iteration,beta0,B0,shape,scale,nu,v,G0,A,amount,customer,process)
paramters=test1[[1]]
variance=test1[[2]]

library(ggplot2)
qplot(1:10000,paramters[,7])
temp=colMeans(paramters[5000:10000,])
write.csv(temp,'temp.csv')

ind=seq((5000*3+3),10000*3,by=3)
mean(variance[ind,3])
##question 1 why don't use hierachical setup for beta
##question 2 there is not interaction between transaction happening and gross margin?