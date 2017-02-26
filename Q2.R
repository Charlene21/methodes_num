#Simulation du processus de Lévy de type diffusion à sauts
SimulModele<-function(lambda, T, mu, delta, sigma, pas,S0,p, lambda1, lambda2,bool){
  
  #on pose gamma tel qu'il est donné en exercice 3
  gamma = lambda * (1-exp(mu+delta*delta/2)) - (sigma*sigma)/2 + 0.01;
  
  #Simulation du processus de Poisson composé
  params = SimulTrajPoissonCompose(lambda, T, mu, delta,p, lambda1, lambda2,bool);
  Tn = params$Tn
  Tn=c(0,Tn,T)
  X_compose = params$PoissonCompose
  
  #Discrétisation
  Tn_recompose = c()
  for (i in 1:(length(Tn)-1)){
    pas_temps = (Tn[i+1]-Tn[i])/pas;
    
    for (j in 0:(pas-1)){
      Tn_recompose = c(Tn_recompose,Tn[i]+j*pas_temps)
    }
    
  }
  
  Tn_recompose = c(Tn_recompose,T)
  
  #On remet les dates dans l'ordre
  Tn_recompose <- sort(Tn_recompose)
  
  X = c()
  j=0;
  
 
  Wt = 0;

  #Simulation du processus de Lévy
  for(i in 1:length(Tn_recompose)){
    
    if(Tn_recompose[i] %in% Tn && i != length(Tn_recompose)){
      j = j+1;
    }

    if(i>1)
      Wt = Wt + sqrt(Tn_recompose[i]-Tn_recompose[i-1]) * rnorm(1, mean=0, sd=1);
    X_t = gamma*(Tn_recompose[i]/T) + sigma*Wt + X_compose[j];
    X = c(X,X_t);
  }
  
  St = c();
  for (i in 1:length(X)){
    St = c(St,S0*exp(X[i]))
  }
  
  
  if (bool == 1){
    plot(Tn_recompose,xlab="Temps", St, main="Modèle de Merton", type='l', panel.first = grid(), lab=c(20,10,12));
  }
  else {
    plot(Tn_recompose,xlab="Temps", St, main="Modèle de Kou", type='l', panel.first = grid(), lab=c(20,10,12)); 
  }
  
}
