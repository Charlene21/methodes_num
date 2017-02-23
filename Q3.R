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
  
  return(exp(-r*T)*gcall(S0*exp((r-(sigma*sigma)/2)*T),sigma*sqrt(T),K));
  
}

seekRoot<-function(N,lambda,T,S0,mu,delta,K,r){
  root = uniroot(function(x) truncCall(N,lambda,T,S0,mu,x,delta,K,r) - BSCall(S0,r,x,T,K), lower = -20, upper = 20, tol=1e-9)$root;
  cat("root : ", root, "\n")
  cat("prix du Call européen Merton: " , truncCall(N,lambda,T,S0,mu,root,delta,K,r), "\n");
  cat("prix du Call européen B&S : " , BSCall(S0,r,root,T,K));
  return(root);
}

traceCourbeVolatilite <- function(N,lambda,T,S0,mu,delta,r){
  volatilite = c();
  strike = c();
  inf = 0.8*S0;
  sup = 1.2*S0;
  for (K in inf:sup){
    cat("K = ", K, "\n");
    strike = c(strike,K);
    volatilite = c(volatilite, seekRoot(N,lambda,T,S0,mu,delta,K,r));
  }
  
  cat("strike : ", strike, "\n");
  cat("vol : ", volatilite, "\n");
  
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