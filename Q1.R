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
  cat("Exp : ", Exp, "\n")
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
SimulTrajPoissonCompose<-function(lambda, T, mu, delta){
  X = c(0);
  somme = 0;
  Tn = SimulTrajPoisson(lambda, T);
  NtFinal = length(Tn);
  
  for (i in 1:NtFinal){
    Y <- rnorm(1, mean=mu, sd=delta);
    cat("Y : ", Y, "\n")
    somme = somme + Y;
    X = c(X,somme);
  }

  cat("X :" , X, "\n")
  plot(c(0,Tn),X, main="compose", type='s'); 
  
  Params <- list(PoissonCompose = X, Tn = Tn)
  return(Params)
}
