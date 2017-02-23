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
}

gcall <- function(x,sigma,K){
  d1 = log(x/K)/sigma;
  d2 = d1 + sigma;
  return (x*exp(sigma*sigma/2)*pnorm(d2) - K*pnorm(d1));
}

truncCall <- function(N,lambda,t, S0, mu, sigma, delta,K, r){
  sum = 0;
  gamma = lambda * (1-exp(mu+delta*delta/2)) - (sigma*sigma)/2 +r;
  for (k in 0:N){
    #cat((lambda*T)^k, "\n")
    terme = (((lambda*t)^k)/factorial(k)) * gcall(S0*exp(gamma*t+k*mu), sqrt(sigma*sigma*t+k*delta*delta),K);
    sum = sum + terme;
  }
  sum = sum * exp(-t*(r+lambda));
}

erreurExacte<-function(S0,K,r,t){
  return (S0 - K*exp(-r*t));
}

#erreurSimul(10,0.5,20, 100,0.4 , 0, 0.25, 0.02,80, 0.01)
erreurSimul<-function(N,lambda,t, S0, mu, sigma, delta,K, r){
  put = truncPut(N,lambda,t, S0, mu, sigma, delta,K, r);
  call = truncCall(N,lambda,t, S0, mu, sigma, delta,K, r);
  cat("put = ", put, "\n");
  cat("call = ", call, "\n");
  return (call - put);
}

comparaisonErreur <- function(N,lambda,t, S0, mu, sigma, delta,K, r){
  erreur_sim = erreurSimul(N,lambda,t, S0, mu, sigma, delta,K, r);
  erreur_exacte = erreurExacte(S0,K,r,t);
  cat("erreur simulée : ", erreur_sim, "\n");
  cat("erreur exacte : ", erreur_exacte, "\n");

}

BSCall <- function(S0, r, sigma, T, K){
  
  return(exp(-r*T)*gcall(S0*exp((r-sigma*sigma/2)*T),sigma*sqrt(T),K));
  
}

seekRoot<-function(N,lambda,T,S0,mu,delta,K,r,sigma){
  root = uniroot(function(x) BSCall(S0,r,x,T,K) - truncCall(N,lambda,T,S0,mu,x,delta,K,r), lower = -20, upper = 20, tol=1e-9)$root;
  cat("root : ", root, "\n")
  cat("prix du Call européen Merton: " , truncCall(N,lambda,T,S0,mu,sigma,delta,K,r), "\n");
  cat("prix du Call européen B&S : " , BSCall(S0,r,sigma,T,K));
}