#########################################################################
#
# Classifier based on HMC and MLE estimation
#
#########################################################################

#install.packages(c("gtools","glmnet","Rfast","arrayhelpers"))
library(gtools)
library(class)
library(glmnet)
library(Rfast)
library(arrayhelpers)

##################
# Functions
##################

#Function to count the occurance of sequences
count.seq<-function(test,my.seq)
{
  n<-length(test)
  m<-length(my.seq)
  count=0
  for(i in 1:(n-m+1))
  {
    if(prod(test[i:(i+m-1)]==my.seq)==1)
    {
      count=count+1
    }
  }
  return(count)
}

#Function to calculate the vector of counting statistics
my.count.stat<-function(test,k)
{
  all.perms<-permutations(2,k,0:1,repeats=T)#for binary sequence only. otherwise change 0:1 
  all.counts<-numeric(nrow(all.perms))
  for(j in 1:nrow(all.perms))
  {
    all.counts[j]<-count.seq(test,all.perms[j,])
  }
  return(all.counts)
}

#Calculation of initial probability
initial.prob<-function(data,k)
{
  return(prop.table(table(data.frame(data[,1:k]))))
}

#Function for simulation of a second order MC
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

#Calculation of counting statistics for training data
freq.calc<-function(data,k)
{
  train.stats<-matrix(0,nrow = nrow(data),ncol=2^(k+1))
  for(i in 1:nrow(train.stats))
  {
    train.stats[i,]<-my.count.stat(data[i,],k+1)
  }
  return(train.stats)
}

#Calculation of the paramerters regrading conditional probability
p.cond_calc<-function(stats)
{
  sum.cond<-colSums(stats)
  
  p.cond<-numeric(length(sum.cond))
  for(i in 1:length(sum.cond))
  {
    if(i%%2==1)
    {
      p.cond[i]<-sum.cond[i]/sum(sum.cond[i:(i+1)])
    }
    else
    {
      p.cond[i]<-sum.cond[i]/sum(sum.cond[(i-1):i])
    }
  }
  return(p.cond)
}

#calculation of log-likelihood
log.lhd<-function(initial.pr,p.cond,seq,k)
{
  if(k==1)
  {
    initial.index<-seq[1:k]+1
  }else
  {
    initial.index<-array2vec(seq[1:k]+1,dim=rep(2,k))
  }
  dum<-log(initial.pr[initial.index])
  
  counts<-my.count.stat(seq,k+1)
  for(i in 1:length(p.cond))
  {
    dum=dum+counts[i]*log(p.cond[i])
  }
  return(dum)
}

#calculation of CV-error for a choices of orders of MC
my.test<-function(train.set.data,label,k0,k1)
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
    train.stats_0<-freq.calc(train.data.cv[train.label.cv==0,],k0)
    train.stats_1<-freq.calc(train.data.cv[train.label.cv==1,],k1)
    
    p.cond_0<-p.cond_calc(train.stats_0)
    p.cond_1<-p.cond_calc(train.stats_1)
    
    initial.prob_0<-initial.prob(train.data.cv[train.label.cv==0,],k0)
    initial.prob_1<-initial.prob(train.data.cv[train.label.cv==1,],k1)
    
    est.label<-numeric(nrow(test.data.cv))
    for(i in 1:nrow(test.data.cv))
    {
      est.label[i]<-(log.lhd(initial.prob_0,p.cond_0,test.data.cv[i,],k0))<(log.lhd(initial.prob_1,p.cond_1,test.data.cv[i,],k1))
    }
    foo<-foo+(sum(test.label.cv!=est.label)/length(test.label.cv))
  }
  return(foo/K)
}

#calculation of test error
test.error.calc<-function(train.set.data,label,test.set.data,label1,k0,k1)
{
  #Final estimation of parameters
  f.train.stats_0<-freq.calc(train.set.data[label==0,],k0)
  f.train.stats_1<-freq.calc(train.set.data[label==1,],k1)
  
  f.p.cond_0<-p.cond_calc(f.train.stats_0)
  f.p.cond_1<-p.cond_calc(f.train.stats_1)
  
  f.initial.prob_0<-initial.prob(train.set.data[label==0,],k0)
  f.initial.prob_1<-initial.prob(train.set.data[label==1,],k1)
  
  est.label<-numeric(nrow(test.set.data))
  for(i in 1:nrow(test.set.data))
  {
    est.label[i]<-(log.lhd(f.initial.prob_0,f.p.cond_0,test.set.data[i,],k0))<(log.lhd(f.initial.prob_1,f.p.cond_1,test.set.data[i,],k1))
  }
  return(sum(label1!=est.label)/length(label1))
}

######################################
# Replication
######################################
R<-100 #No. of replications
results<-matrix(0,nrow = R,ncol=4)
r<-1
count<-1

while(count<=R)
{
  start.time<-Sys.time()
  ####################
  # Data Simulation
  ####################
  #Generation of Training Data
  set.seed(r) 
  k0=2;k1=2
  
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
  
  ##################
  #Calculations
  ##################
  #Calculation of cross-valudation errors for various choices of orders of MC used
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
  
  order<-numeric(2)
  if(nrow(est)==1)
  {
    order<-est
  }
  else
  {
    order<-est[which(rowSums(est)==min(rowSums(est))),]
  }
  
  print(order)
  dum<-test.error.calc(train.set.data,label,test.set.data,label1,order[1],order[2])
  
  end.time<-Sys.time()
  
  results[count,]<-c(order[1],order[2],dum,end.time-start.time)
  count=count+1
  r=r+1
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

write.table(results,file="MLE.txt",row.names = F, col.names = F)
#write.csv(results,file="MLE.csv",row.names = F, col.names = F)