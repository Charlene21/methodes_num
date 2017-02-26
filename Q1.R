#question a
SimulTrajPoisson<-function(lambda, T)
{
  somme =0;
  
  NtFinal = 0;
  Tau = c()
  Tn_while = c();
  Tn = c();
  Nt = c();
  
  while (somme < T){
    
    #génération du vecteur de variables de loi exponentielle
    Exp = rexp(n=1, rate = lambda);  #rate = intensité
    Tau = c(Tau, Exp);
    
    
    somme = somme + Exp;
    
    Tn_while = c(Tn_while , somme);
    
  }
  
  nr = length(Tn_while)

  if(length(Tn_while) > 1){
    
    for (i in 1:(nr-1)){
      Tn = c(Tn,Tn_while[i])
    }
    
    
    #Génération des variables de poisson
    for (j in 1:(nr-1)){
      if (Tn[j] <= T) {
        NtFinal = NtFinal + 1;
        Nt[j] = NtFinal;
      }
    }
    
    plot(c(0,Tn, T),c(0:length(Tn), length(Tn)), xlab = "Temps", ylab="NT",  type='s', main="Processus de Poisson"); 
  }
  else {
    cat("Il n'y a pas de temps de saut \n ");
  }
  
  return(Tn)
}



#question b
SimulTrajPoissonCompose<-function(lambda, T, mu, delta,p, lambda1, lambda2,bool){
  X = c(0);
  somme = 0;
  Tn = SimulTrajPoisson(lambda, T);
  
  if (length(Tn) != 0){
  NtFinal = length(Tn);
  
  for (i in 1:NtFinal){
    if (bool == 1){
      Y <- rnorm(1, mean=mu, sd=delta);
    }
    
    else{
      Y <- generateKou(p, lambda1, lambda2);
    }
    
    somme = somme + Y;
    X = c(X,somme);
  }
  

  plot(c(0,Tn,T),c(X,X[length(X)]), xlab = "Temps", ylab= "Xt", main="Processus de Poisson composé", type='s'); 

  }
  
  
  Params <- list(PoissonCompose = X, Tn = Tn)
  return(Params)
}

generateKou <- function(p, lambda1, lambda2){
  U1 = runif(1,0,1)
  U2 = runif(1,0,1)
  if(U1 <= p){
    x = log(U2);
  } else {
    x= -log(U2);
  }
  
  if(x>0){
    Y <- p*lambda1*exp(-lambda1*x);
  } else {
    Y <- (1-p)*lambda2*exp(-lambda2*abs(x));
  }
  
  return(Y);
}
