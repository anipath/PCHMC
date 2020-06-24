#############################################################################################
#
#Classifer based on Mixture Models of Higher Order Chains
#
#############################################################################################

#install.packages(c("markovchain","arrayhelpers"))
library(markovchain)
library(arrayhelpers)

##################
# Functions
##################

#Function for simulation of a second order MC
MC_1<-function(n,M,initial)
{
  x<-numeric(n)
  x[1]<-sample(x=c(0,1),1,replace = T,prob=M[initial+1,])
  for(i in 2:n)
  {
    x[i]<-sample(x=c(0,1),1,replace = T,prob=M[x[i-1]+1,])
  }
  return(x)
}

MC_2<-function(n,M,initial)
{
  x<-numeric(n)
  x[1]<-sample(x=c(1,2),1,replace = T,prob=M[initial[1]+1,initial[2]+1,])
  x[2]<-sample(x=c(1,2),1,replace = T,prob=M[initial[2]+1,x[1],])
  for(i in 3:n)
  {
    x[i]<-sample(x=c(1,2),1,replace = T,prob=M[x[i-2],x[i-1],])
  }
  x[x==1]=0;x[x==2]=1;
  return(x)
}

#Functions to calculation the probabilities
prob.cond.mc<-function(model,x,given)
{
  lambda<-model$avg.lambda
  n<-length(lambda)
  p<-numeric(n)
  for(i in 1:n)
  {
    p[i]<-model$avg.Q[[i]][x+1,given[n-i+1]+1]
  }
  return(as.numeric(crossprod(lambda,p)))
}

initial.prob<-function(data,k)
{
  return(prop.table(table(data.frame(data[,1:k]))))
}

log.prob.mc.seq<-function(data,model,x)
{
  k<-length(model$avg.lambda)
  m<-length(unique(data[1,]))
  n<-length(x)
  if(k==1)
  {
    initial.index<-x[1:k]+1
  }else
  {
    initial.index<-array2vec(x[1:k]+1,dim=rep(m,k))
  }
  dum<-log(initial.prob(data,k)[initial.index])
  for(i in 1:(n-k))
  {
    dum=dum+log(prob.cond.mc(model,x[i+k],x[i:(i+k-1)]))
  }
  return(dum)
}

#Estimation of Model Parameters
my.fit.hmc<-function(data,k)
{
  m<-length(unique(data[1,]))
  store<-list(avg.lambda=numeric(k),avg.Q=rep(list(matrix(0,m,m)),k),avg.X=numeric(m))
  n<-nrow(data)
  for(i in 1:n)
  {
    model<-fitHigherOrder(data[i,],order = k)
    
    store$avg.lambda=store$avg.lambda+model$lambda
    for(j in 1:k)
    {
      store$avg.Q[[j]]=store$avg.Q[[j]]+model$Q[[j]]
    }
    store$avg.X=store$avg.X+model$X
  }
  store$avg.lambda=store$avg.lambda/n
  for(j in 1:k)
  {
    store$avg.Q[[j]]=store$avg.Q[[j]]/n
  }
  store$avg.X=store$avg.X/n
  return(store)
}

my.test<-function(train.set.data,label,order.1,order.2)
{
  K <- 2
  n<-nrow(train.set.data)
  permutation <- sample(1:n, replace = FALSE)
  test.index <- split(permutation, rep(1:K, length = n, each = n/K))
  
  foo <- 0
  for(k in 1:K)
  {
    train.data.cv <- train.set.data[-test.index[[k]], ]
    train.label.cv <- label[-test.index[[k]]]
    test.data.cv <- train.set.data[test.index[[k]], ]
    test.label.cv <- label[test.index[[k]]]
    
    #Fitting of models
    model.0<-my.fit.hmc(train.data.cv[train.label.cv==0,],order.1)
    model.1<-my.fit.hmc(train.data.cv[train.label.cv==1,],order.2)
    
    est.label<-numeric(nrow(test.data.cv))
    for(i in 1:nrow(test.data.cv))
    {
      est.label[i]<-(log.prob.mc.seq(train.data.cv[train.label.cv==0,],model.0,test.data.cv[i,])<log.prob.mc.seq(train.data.cv[train.label.cv==1,],model.1,test.data.cv[i,]))
    }
    foo<-foo+(sum(test.label.cv!=est.label)/length(test.label.cv))
  }
  return(foo/K)
}

#calculation of test error
test.error.calc<-function(train.set.data,label,test.set.data,label1,est)
{
  #Fitting of final models
  model.0<-my.fit.hmc(train.set.data[label==0,],est[1])
  model.1<-my.fit.hmc(train.set.data[label==1,],est[2])
  
  est.label<-numeric(nrow(test.set.data))
  for(i in 1:nrow(test.set.data))
  {
    est.label[i]<-(log.prob.mc.seq(train.set.data[label==0,],model.0,test.set.data[i,])<log.prob.mc.seq(train.set.data[label==1,],model.1,test.set.data[i,]))
  }
  
  return(sum(label1!=est.label)/length(label1))
}


###################
# Data Simulation
###################
R<-100 #No. of replications
results<-matrix(0,nrow = R,ncol=4)
count<-1
r<-1

while(count<=R)
{
  start.time<-Sys.time()
  ####################
  # Data Simulation
  ####################
  #Generation of Training Data
  set.seed(r) 
  k0=2;k1=2;
  
  M0=array(c(0.4,0.6,0.6,0.4,0.6,0.4,0.4,0.6),dim = c(2,2,2))
  M0_1<-array(c(0.4,0.6,0.6,0.4),dim = c(2,2))
  train0<-matrix(nrow = 50,ncol = 100)
  for(i in 1:50)
  {
    if(k0==2)
    {
      train0[i,]<-MC_2(100,M0,c(0,1))
    }
    else
    {
      train0[i,]<-MC_1(100,M0_1,0)
    }
  }
  
  train1<-matrix(nrow = 50,ncol = 100)
  M1=array(c(0.6,0.4,0.4,0.6,0.4,0.6,0.6,0.4),dim = c(2,2,2))
  M1_1<-array(c(0.6,0.4,0.4,0.6),dim = c(2,2))
  for(i in 1:50)
  {
    if(k1==2)
    {
      train1[i,]<-MC_2(100,M1,c(0,1))
    }
    else
    {
      train1[i,]<-MC_1(100,M1_1,0)
    }
  }
  
  label<-c(rep(0,50),rep(1,50))
  train.set.data<-rbind(train0,train1)
  
  #Generation of Test Data
  test0<-matrix(nrow = 150,ncol = 100)
  for(i in 1:150)
  {
    if(k0==2)
    {
      test0[i,]<-MC_2(100,M0,c(0,1))
    }
    else
    {
      test0[i,]<-MC_1(100,M0_1,0)
    }
  }
  
  test1<-matrix(nrow = 150,ncol = 100)
  for(i in 1:150)
  {
    if(k1==2)
    {
      test1[i,]<-MC_2(100,M1,c(0,1))
    }
    else
    {
      test1[i,]<-MC_1(100,M1_1,0)
    }
  }
  
  label1<-c(rep(0,150),rep(1,150))
  test.set.data<-rbind(test0,test1)
  
  ################
  # Calculations
  ################
  
  miss.result<-matrix(0,3,3)
  for(i in 1:3)
  {
    for(j in 1:3)
    {
      dum<-0
      for(k in 1:10)
      {
        dum<-dum+my.test(train.set.data,label,i,j)
      }
      miss.result[i,j]<-dum/10
      print(miss.result[i,j])
    }
  }
  
  if(sum(is.na(miss.result)!=0))
  {
    r=r+1
    next
  }
  
  colnames(miss.result)<-1:3;rownames(miss.result)<-1:3
  print(round(miss.result,4))
  
  est<-which(miss.result==min(miss.result),arr.ind = TRUE)
  print(est)
  
  order<-numeric(2)
  if(nrow(est)==1)
  {
    order<-est
  }
  else
  {
    order<-est[which(rowSums(est)==min(rowSums(est))),]
  }
  
  dum<-test.error.calc(train.set.data,label,test.set.data,label1,order)
  print(dum)
  
  end.time<-Sys.time()
  
  results[count,]<-c(order[1],order[2],dum,end.time-start.time)
  r=r+1
  count=count+1
}

results

table(data.frame(results[,1:2]))/sum(table(data.frame(results[,1:2])))

mean(results[,3])
sd(results[,3])

mean(results[,4])
sd(results[,4])

output<-list(
  results,
  table(data.frame(results[,1:2]))/sum(table(data.frame(results[,1:2]))),
  mean(results[,3]),sd(results[,3]),
  mean(results[,4]),sd(results[,4])
)

write.table(results,file="HMC.txt",row.names = F, col.names = F)
#write.csv(results,file="HMC.csv",row.names = F, col.names = F)