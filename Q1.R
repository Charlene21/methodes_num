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
  #cat("Tn_while : ", Tn_while, "\n")
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

  cat("Nt : ", Nt, "\n")


plot(c(0,Tn),0:length(Tn), type='s', main="poisson"); 
  
return(Tn)
}



#question b
SimulTrajPoissonCompose<-function(lambda, T, mu, delta,p, lambda1, lambda2,bool){
  X = c(0);
  somme = 0;
  Tn = SimulTrajPoisson(lambda, T);
  NtFinal = length(Tn);
  
  for (i in 1:NtFinal){
    if (bool == 1){
      Y <- rnorm(1, mean=mu, sd=delta);
    }
      
    else{
      Y <- generateKou(p, lambda1, lambda2);
    }
      
    #cat("Y : ", Y, "\n")
    somme = somme + Y;
    X = c(X,somme);
  }

  #cat("X :" , X, "\n")
  plot(c(0,Tn),X, main="compose", type='s'); 
  
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
