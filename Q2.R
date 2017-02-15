#Question a
source('/user/9/.base/paviotch/home/3A/scilab/Q1.R')

ModeleMerton<-function(lambda, T, mu, delta, gamma, sigma, pas){
  params = SimulTrajPoissonCompose(lambda, T, mu, delta);
  Tn = params$Tn;
  
  discretisation = seq(0,T,pas);

  #cat("discretisation : " , discretisation)
  discretisation = c(discretisation,Tn)

  #cat("discretisation : " , discretisation, "\n")
  discretisation <- sort(discretisation)
  #cat("discretisation : " , discretisation, "\n")

  j = 1;
  X = c();
  #cat("length : " , length(discretisation), " / " , discretisation[2], "\n")
  
  for(i in 1:length(discretisation)){
   # cat("i = " , discretisation[i] )
    params = SimulTrajPoissonCompose(lambda, T, mu, delta);
    PoissonCompose = params$PoissonCompose
    
    W <- rnorm(1, mean=0, sd=discretisation[i]);
    X[j] = gamma*discretisation[i] + sigma*W + PoissonCompose;
    j = j+1;
  }
  
  cat("X : " , X, "\n")
  return(X)
}