#### GIBBS SAMPLING FOR THE UNINITIATED
require(MCMCpack)


#### function to simulate truth
simu=function(N=100, V=5, gamma.pi.1=2, gamma.pi.0=8, gamma.theta=rep(1,V), n.words=c(100,200)){
  
  #gamma.pi.1=2
  #gamma.pi.0=8
  
  
  #### generate true values:
  ## Label of documents: 
  # pi, prob in the bernoulli trail
  pi=rbeta(N, gamma.pi.1, gamma.pi.0)
  # label L
  L=rbinom(N, 1, pi)
  
  ## Bag of words in each documents:
  # prabability of words, for L=0 and L=1
  theta=rdirichlet(2, gamma.theta)
  
  # number of words in documents; suppose 100-200 words
  R=sample(n.words[1]:n.words[2], size=N, replace = T)
  
  # Bag of words
  W=NULL
  for(i in 1:N){
    if(L[i]==1){
      theta.i=theta[2,]
    }else{
      theta.i=theta[1,]
    }
    W=cbind(W,rmultinom(n=1, size=R[i], prob=theta.i))
  }
  
  res=list(L=L, W=W, R=R, pi=pi, theta=theta, gamma.pi.1=gamma.pi.1, gamma.pi.0=gamma.pi.0)
  res
}

#### function for gibbs sampler
gibbs=function(N,V,W, iter=10){ ## N: number of 
  ### (1) init sampler
  
  #current label
  L=rbinom(N, 1, 0.5) #just random guess
  
  #current value of theta_1 and theta_0
  theta.0=rdirichlet(1,rep(1,V))
  theta.1=rdirichlet(1,rep(1,V))
  
  # number of docs labeled as 1/0
  c.1=sum(L)
  c.0=N-c.1
  
  # word counts for all docs labeled as 1/0
  w.counts.1=rowSums(W[,which(L==1)])
  w.counts.0=rowSums(W[,which(L==0)])
  
  # non-info prior
  gamma.pi.1=1
  gamma.pi.0=1
  
  ### (2) start sampling
  
  ## keep record of paras
  L.all=L
  theta.0.all=theta.0
  theta.1.all=theta.1
  
  for( i in 1:iter){
    for(j in 1:N){
      w.count.j=W[,j]
      
      #change C_1/C_0 and w.counts
      if(L[j]==1){
        w.counts.1 = w.counts.1 - w.count.j
        c.1 = c.1 - 1
      }else{
        w.counts.0 = w.counts.0 - w.count.j
        c.0 = c.0 - 1
      }
      
      #assign a new label
      log.val.0= log(c.0 + gamma.pi.0 - 1) - log(N + gamma.pi.1 + gamma.pi.0 - 1) + sum(log(theta.0) * w.count.j)
      log.val.1= log(c.1 + gamma.pi.1 - 1) - log(N + gamma.pi.1 + gamma.pi.0 - 1) + sum(log(theta.1) * w.count.j)
      rr=exp(log.val.1 - log.val.0)
      if(is.infinite(rr)){
        pp=1
      }else{
        pp=rr/(1+rr)
      }
      L.j.next=rbinom(1,1, pp)
      L[j] = L.j.next
      
      #update counts
      if(L.j.next==1){
        c.1 = c.1 + 1
        w.counts.1 = w.counts.1 + w.count.j
      }else{
        c.0 = c.0 + 1
        w.counts.0 = w.counts.0 + w.count.j
      }
      
    }
    
    theta.0=rdirichlet(1, w.counts.0 + gamma.pi.0)#vector of total word counts, including pseudocounts
    theta.1=rdirichlet(1, w.counts.1 + gamma.pi.1)
    
    #record
    L.all=rbind(L.all, L)
    theta.0.all=rbind(theta.0.all, theta.0)
    theta.1.all=rbind(theta.1.all, theta.1)
    
    if(i %% 10 ==0){
      message(paste0("iter=",i))
    }
  }
  
  per.iter=list(theta.0=theta.0.all, theta.1=theta.1.all, L=L.all)
  theta=rbind(theta.0, theta.1)
  res=list(theta=theta, per.iter=per.iter, iter=iter)
  res
  
}










