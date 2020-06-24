############################################################################################
#
# LASSO-Logistic Classifier based on Counting based statistics
# 
############################################################################################

#install.packages(c("gtools","glmnet","Rfast"))
library(gtools)
library(class)
library(glmnet)
library(Rfast)


############################
#Functions
############################

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
my.test<-function(train.set.data,label,k)
{
  train.stats<-matrix(0,nrow = nrow(train.set.data),ncol=2^k)
  for(i in 1:nrow(train.stats))
  {
    train.stats[i,]<-my.count.stat(train.set.data[i,],k)
  }
  
  K <-2
  n<-nrow(train.stats)
  permutation <- sample(1:n, replace = FALSE)
  test.index <- split(permutation, rep(1:K, length = n, each = n/K))
  
  foo <- 0
  for(k in 1:K)
  {
    train.data.cv <- train.stats[-test.index[[k]], ]
    train.label.cv <- label[-test.index[[k]]]
    test.data.cv <- train.stats[test.index[[k]], ]
    test.label.cv <- label[test.index[[k]]]
    
    
    #Fitting of GLMNET model
    glm.model = glmnet(as.matrix(train.data.cv), train.label.cv, family = "binomial")
    l.range<-seq(0,0.35,by=0.01)
    prd<-predict(glm.model, newx = as.matrix(test.data.cv),type="response",s=l.range)
    final.prd<-round(prd)
    #Calculation of Misclassification
    dum<-numeric(ncol(final.prd))
    for(i in 1:ncol(final.prd))
    {
      dum[i]<-sum(final.prd[,i]!=test.label.cv)/length(test.label.cv)
    }
    foo<-foo+dum
  }
  if(sum(is.na(foo))!=0)
  {
    next
  }
  return(foo/K)
}

######################################
# Replication
######################################

R<-100 #No. of replications
results<-matrix(0,nrow = R,ncol=4)
for(r in 1:R)
{
  start.time<-Sys.time()
  ######################################
  # Data Simulation
  ######################################
  
  set.seed(r)
  #Generation of Training Data
  k=2
  M0=array(c(0.4,0.6,0.6,0.4,0.6,0.4,0.4,0.6),dim = c(2,2,2))
  M0_1<-array(c(0.4,0.6,0.6,0.4),dim = c(2,2))
  train0<-matrix(nrow = 50,ncol = 100)
  for(i in 1:50)
  {
    if(k==2)
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
    if(k==2)
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
    test0[i,]<-MC_2(100,M0,c(0,1))
  }
  test1<-matrix(nrow = 150,ncol = 100)
  for(i in 1:150)
  {
    test1[i,]<-MC_2(100,M1,c(0,1))
  }
  label1<-c(rep(0,150),rep(1,150))
  test.set.data<-rbind(test0,test1)
  
  
  #############################
  #Calculations
  #############################
  
  #Calculation and Plot of CV-errors 
  l.range<-seq(0,0.35,by=0.01)
  miss.error<-matrix(0,nrow = 5,ncol = length(l.range))
  
  dum<-numeric(length(l.range))
  for(i in 1:10)
  {
    dum= dum+my.test(train.set.data,label,1)
  }
  miss.error[1,]<-dum/10
  
  for(k in 2:5)
  {
    dum<-numeric(length(l.range))
    for(i in 1:10)
    {
      dum=dum+my.test(train.set.data,label,k)
    }
    miss.error[k,]<-dum/10
  }
  
  #Minimization of CV error over lambda
  miss.error.mins<-rowMins(miss.error,value = T)
  
  print(miss.error.mins)
  
  miss.error.min.lambda<-l.range[rowMins(miss.error)]
  
  #minimization of CV error over k
  est<-list(k=min(which(miss.error.mins==min(miss.error.mins))),lambda=miss.error.min.lambda[k])
  print(est)
  
  #Calculation of counting statistics for test data
  k=est$k
  train.stats<-matrix(0,nrow = nrow(train.set.data),ncol=2^k)
  for(i in 1:nrow(train.stats))
  {
    train.stats[i,]<-my.count.stat(train.set.data[i,],k)
  }
  
  test.stats<-matrix(nrow = nrow(test.set.data),ncol=2^k)
  for(i in 1:nrow(test.stats))
  {
    test.stats[i,]<-my.count.stat(test.set.data[i,],k)
  }
  
  #Final model
  glm.model = glmnet(as.matrix(train.stats), label, family = "binomial")
  l.range<-seq(0,0.35,by=0.01)
  prd<-predict(glm.model, newx = as.matrix(test.stats),type="response",s=est$lambda)
  final.prd<-round(prd)
  
  test.result<-numeric(ncol(final.prd))
  for(k in 1:ncol(final.prd))
  {
    test.result[k]<-sum(final.prd[,k]!=label1)/length(label1)
  }
  
  #Test error
  print(test.result)
  end.time<-Sys.time()
  
  results[r,]<-c(est$k,est$lambda,test.result,end.time-start.time)
}

results

table(results[,1])/R

mean(results[,3])
sd(results[,3])

mean(results[,4])
sd(results[,4])

output<-list(
  full_results=results,
  EPD=table(data.frame(results[,1:2]))/sum(table(data.frame(results[,1:2]))),
  avg_error=mean(results[,3]),sd_error=sd(results[,3]),
  avg_time=mean(results[,4]),sd_time=sd(results[,4])
)

write.table(results,file="count.txt",row.names = F, col.names = F)
#write.csv(results,file="count.csv",row.names = F, col.names = F)