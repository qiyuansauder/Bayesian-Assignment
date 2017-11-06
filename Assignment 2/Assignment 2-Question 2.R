######################################
###COMM 590  Assignment 1 Question 2
###Qiyuan Wang
######################################
library(mvtnorm)

#Metropolis-Hastings Algorithm for Multinomial Logit Model
#in this model we assume number of alternatives in each choice is the same
#and there is not outside option

#choice probability
PChoice=function(X,theta,Choice){
  #number of choices made in the data
  noption=nrow(X)/length(Choice)
  nchoice=length(Choice)
  #utility of each option
  Utility=exp(X%*%theta)
  #reshape into nchoice*noption matrix
  temp=matrix(Utility,ncol=noption,byrow=T)
  #calculate the choice probability
    #the index for choose option
    ind=cbind(1:nchoice,Choice)
    #do log transformation here because the product of 1000 probability will be 0
    choiceprob=sum(log(temp[ind]/(rowSums(temp))))
  
  return(choiceprob)
}

#prior probability
PPrior=function(theta,mean,variance){
  PPrior=log(dmvnorm(theta,mean=mean,sigma=variance))
  return(PPrior)
}

#posterior probability
targetRelProb=function(X,theta,Choice,mean,variance){
  PPosterior=PChoice(X,theta,Choice)+PPrior(theta,mean,variance)
  return(PPosterior)
}


#the main function
MultiNomial=function(X,Choice,N,mean,variance){
  #X is independent variable
  #Y is dependent variable, namely which option is chosen
  #mean and variance is the parameter for normal posterior

  #number of parameters
  nparams=ncol(X)
  
  #construct matrix to store iteration results
  draw=matrix(0,nrow=N,ncol=nparams)
  #set the starting value
  draw[1,]=1
  
  thetaold=draw[1,]
  #the iteration for HM algorithm
  for(i in 2:N){
    #loop over the parameters
    for(j in 1:nparams){
      #using random walk propose
      thetaprime=rnorm(1,mean=0,sd=2)
      #proposed theta
      thetanew=thetaold
      thetanew[j]=thetanew[j]+thetaprime
      
      #calculate the probability of accepting
      PAccept=min(c(1,exp(targetRelProb(X,thetanew,Choice,mean,variance)-
                      targetRelProb(X,thetaold,Choice,mean,variance))))
      
      if(runif(1)<=PAccept){
        thetaold=thetanew
      }
      else{
        thetanew=thetaold
      }
    }
    draw[i,]=thetanew
    if(i<=100){
      cat("Iteration",i,": ",draw[i,], "\n")
    }else{
      #every 100 iteration print the results
      if(i%%100==0){
        cat("Iteration",i,": ",draw[i,], "\n")
      }
    }
  }
  return(draw)
}


###test the model
#generate data
#In the dataset, there are 3 IV (including intercept).
#For each choice, there are 3 options
theta=c(-1,0.5,1)
X=matrix(rnorm(6000),ncol=2)

#calculate utility and choice probability
temputility=exp(X%*%theta)
temputility=matrix(temputility,ncol=3,byrow=T)
tempPchoice=temputility/rowSums(temputility)

#make the choice
Choice=numeric(1000)
tempPchoice[,2]=tempPchoice[,2]+tempPchoice[,1]
tempP=runif(1000)
for (i in 1:1000){
  if(tempP[i]<=tempPchoice[i,1]){
    Choice[i]=1
  }
  else if(tempP[i]<=tempPchoice[i,2]){
    Choice[i]=2
  }
  else{
    Choice[i]=3
  }
}

#run the code
mean=c(1,1)
variance=diag(c(1,1))
N=10000
test=MultiNomial(X,Choice,N,mean,variance)

