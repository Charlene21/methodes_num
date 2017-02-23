#Question a
#SimulModele(0.5,20,0,0.02,0.02,50,100,0.3,0.2,0.5,0)
SimulModele<-function(lambda, T, mu, delta, sigma, pas,S0,p, lambda1, lambda2,bool){
  gamma = 0; #lambda * (1-exp(mu+delta*delta/2)) - (sigma*sigma)/2 + 0.01;
  cat("gamma : ", gamma);
  params = SimulTrajPoissonCompose(lambda, T, mu, delta,p, lambda1, lambda2,bool);
  Tn = params$Tn
  Tn=c(0,Tn,T)
  X_compose = params$PoissonCompose
  #cat("length : ", length(Tn), "\n")
  cat("Tn : " , Tn, "\n")
  #cat("X_compose : ", X_compose, "\n")
  
  Tn_recompose = c()
  for (i in 1:(length(Tn)-1)){
    #cat("Tn[i+1]-Tn[i] : ", Tn[i+1], "-", Tn[i], "\n")
    pas_temps = (Tn[i+1]-Tn[i])/pas;
    #cat("pas temps : ", pas_temps, "\n")
    for (j in 0:(pas-1)){
      Tn_recompose = c(Tn_recompose,Tn[i]+j*pas_temps)
    }
    
  }
  
  Tn_recompose = c(Tn_recompose,T)
  Tn_recompose <- sort(Tn_recompose)
  #cat("Tn_recompose : ", Tn_recompose, "\n");
  
  X = c()
  j=0;
  
 
  cat("length : " , length(Tn_recompose), "\n");
  Wt = 0;
  W = c();
  for(i in 1:length(Tn_recompose)){
    #cat("--->", Tn_recompose[i],"\n")
    #cat("test : ", Tn_recompose[i] %in% Tn , "\n")
    if(Tn_recompose[i] %in% Tn && i != length(Tn_recompose)){
      j = j+1;
    }
    #cat("j : ", j, "\n");
    if(i>1)
      Wt = Wt + sqrt(Tn_recompose[i]-Tn_recompose[i-1]) * rnorm(1, mean=0, sd=1);
    
    W=c(W,Wt);
    #cat("Wt : ", Wt, "\n tn_rec : ", Tn_recompose, "\n x_comp : ", X_compose, "\n");
    
    X_t = gamma*(Tn_recompose[i]/T) + sigma*Wt + X_compose[j];
    #cat("X_t : ", X_t, "\n")
   # cat("X_compose[j] : i = ",i," : ", X_compose[j], "\n");
    #cat("Xt = ", X_t, "\n")
    X = c(X,X_t);
  }
  

 # cat("X : ", X, "\n");
  
  St = c()
  for (i in 1:length(X)){
    St = c(St,S0*exp(X[i]))
  }
  
  cat("W : ", W, "\n");
  if (bool == 1){
    plot(Tn_recompose,xlab="Temps", St, main="Modèle de Merton", type='l', panel.first = grid(), lab=c(20,10,12));
  }
  else {
    plot(Tn_recompose,xlab="Temps", St, main="Modèle de Kou", type='l', panel.first = grid(), lab=c(20,10,12)); 
  }
 
  plot(W,xlab="Temps", St, main="Modèle de Merton", type='l', panel.first = grid(), lab=c(20,10,12));
  
  
  cat("Graphique tracé")
}
