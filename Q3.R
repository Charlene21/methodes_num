gput <- function(x,sigma,K){
  d1 = log(K/x)/sigma;
  d2 = d1 - sigma;
  return (K*pnorm(d1) - x*exp(sigma*sigma/2)*pnorm(d2));
}

truncPut <- function(N,lambda,t, S0, mu, sigma, delta,K, r){
  sum = 0;
  gamma = lambda * (1-exp(mu+delta*delta/2)) - (sigma*sigma)/2 +r;
  for (k in 0:N){
    #cat("x = ",S0*exp(gamma*t+k*mu), "\n");
    terme = (((lambda*t)^k)/factorial(k)) * gput(S0*exp(gamma*t+k*mu), sqrt(sigma*sigma*t+k*delta*delta),K);
    sum = sum + terme;
  }
  sum = sum * exp(-t*(r+lambda));
  return(sum);
}

gcall <- function(x,sigma,K){
  d1 = (log(x/K))/sigma;
  d2 = d1 + sigma;
  return (x*exp((sigma*sigma)/2)*pnorm(d2) - K*pnorm(d1));
}

truncCall <- function(N,lambda,T, S0, mu, sigma, delta,K, r){
  sum = 0;
  gamma = lambda * (1-exp(mu+(delta*delta)/2)) - (sigma*sigma)/2 +r;
  for (k in 0:N){
    #cat((lambda*T)^k, "\n")
    terme = (((lambda*T)^k)/factorial(k)) * gcall(S0*exp(gamma*T+k*mu), sqrt(sigma*sigma*T+k*delta*delta),K);
    sum = sum + terme;
  }
  sum = sum * exp(-T*(r+lambda));
  
  return(sum);
}

erreurExacte<-function(S0,K,r,t){
  return (S0 - K*exp(-r*t));
}

erreurSimul<-function(N,lambda,t, S0, mu, sigma, delta,K, r){
  put = truncPut(N,lambda,t, S0, mu, sigma, delta,K, r);
  call = truncCall(N,lambda,t, S0, mu, sigma, delta,K, r);
  cat("N = ", N, " call = ", call, "\nPut : ", put, "\n");
  return (call - put);
}

comparaisonErreur <- function(lambda,T, S0, mu, sigma, delta,K, r, tolerance = 0.00001){
  N = 0;
  erreur_sim = erreurSimul(N,lambda,T, S0, mu, sigma, delta,K, r);
  erreur_exacte = erreurExacte(S0,K,r,T);
  diff = erreur_sim - erreur_exacte;
  
  while(abs(diff) > tolerance){
    N = N+1;
    erreur_sim = erreurSimul(N,lambda,T, S0, mu, sigma, delta,K, r);
    erreur_exacte = erreurExacte(S0,K,r,T);
    diff = erreur_sim - erreur_exacte;
    
  }
  
  call = truncCall(N,lambda,T,S0,mu,sigma,delta,K,r);
  put = truncPut(N,lambda,T,S0,mu,sigma,delta,K,r);
  cat("erreur simulée : ", erreur_sim, "\n");
  cat("erreur exacte : ", erreur_exacte, "\n");
  
  return(N);

}

BSCall <- function(S0, r, sigma, T, K){
  
  return(exp(-r*T)*gcall(S0*exp((r-(sigma*sigma)/2)*T),sigma*sqrt(T),K));
  
}

seekRoot<-function(lambda,T,S0,mu,delta,K,r,sigma, tolerance = 0.00001){
  N = comparaisonErreur(lambda,T, S0, mu, sigma, delta,K, r, tolerance);
  root = uniroot(function(x) truncCall(N,lambda,T,S0,mu,sigma,delta,K,r) - BSCall(S0,r,x,T,K), lower = -2, upper = 2, tol=1e-9)$root;
  cat("root : ", root, "\n")
  cat("prix du Call européen Merton: " , truncCall(N,lambda,T,S0,mu,root,delta,K,r), "\n");
  cat("prix du Call européen B&S : " , BSCall(S0,r,root,T,K));
  return(root);
}

traceCourbeVolatilite <- function(lambda,T,S0,mu,delta,r, sigma){
  volatilite = c();
  strike = c();
  inf = 0.8*S0;
  sup = 1.2*S0;
  for (K in inf:sup){
    cat("K = ", K, "\n");
    strike = c(strike,K);
    volatilite = c(volatilite, seekRoot(lambda,T,S0,mu,delta,K,r,sigma));
  }
  
  plot(x=strike,y=volatilite,main="Volatilité implicite en fonction du strike", type="l");
}

traceCourbe <- function(N,lambda,T,K,mu,sigma,delta,r){
  inf = 10;
  sup = 150;
  call = c();
  put = c();
  cours = c();
  bsCall = c();
  for (S0 in inf:sup){
    bsCall = c(bsCall,BSCall(S0, r, sigma, T, K));
    cours = c(cours, S0);
    call = c(call,truncCall(N,lambda,T,S0,mu,sigma,delta,K,r));
    put = c(put,truncPut(N,lambda,T,S0,mu,sigma,delta,K,r));
  }
  plot(x=cours,y=call,main="Prix des call et put en fonction de S0", type="l", col='green');
  lines(put, col='red');
  lines(bsCall,col='blue');
  legend(20, 90, legend=c("call Merton", "put Merton", "call BS"),
         col=c("green", "red", "blue"),lty=1, cex=0.8)
}

grapheDiff<-function(N,lambda,T,S0,mu,delta,K,r){
  diff = c();
  volatilite = c();
  for (x in seq(0,10,0.1)){
    volatilite=c(volatilite,x);
    diff = c(diff,truncCall(N,lambda,T,S0,mu,x,delta,K,r) - BSCall(S0,r,x,T,K));
  }
  cat("diff : ", diff, "\n");
  plot(x=volatilite,y=diff,main="Ecart de l'erreur", type="l");
}