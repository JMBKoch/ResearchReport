// 2 Factor Model, Cross-Loadings regularized with Small Variance Normal Prior

data{
  int<lower=0> N; // Sample Size
  int<lower=1> P; // Number of Outcomes/ items
  int<lower=1> Q; // Number of Factor
  matrix[N, P] Y; // outcome matrix
  real<lower=0> sigma; // hyperparameter prior cross-loadings
}

parameters{
  vector<lower=0>[P] theta;
  vector[P] lambdaMain;
  vector[P] lambdaCross;
  real<lower=-1,upper=1> factCor;
}

transformed parameters{
  matrix[P, Q] LambdaUnc;
  matrix[P, P] Theta;
  corr_matrix[Q] Psi;
  matrix[P, P] Sigma;
  vector[P] mu;

  Theta = diag_matrix(theta);
  mu = rep_vector(0, P);
  
  // make loading matrix manually; TODO: automate
  LambdaUnc[1:3, 1] = lambdaMain[1:3];
  LambdaUnc[4:6, 2] = lambdaMain[4:6];
  
  LambdaUnc[4:6, 1] = lambdaCross[1:3];
  LambdaUnc[1:3, 2] = lambdaCross[4:6];

  
  // make Psi manually; TODO: automate
  Psi[1, 1] = 1;
  Psi[2, 2] = 1;
  Psi[1, 2] = factCor;
  Psi[2, 1] = factCor;

  Sigma = LambdaUnc*Psi*LambdaUnc' + Theta; 
}

model{
 //priors
 lambdaMain ~ normal(0, 5);
 lambdaCross ~ normal(0, sigma);
 theta ~ cauchy(0, 5);
 // default uniform prior on factCor
  
 //model
 for(i in 1:N)
   Y[i,] ~ multi_normal(mu, Sigma);
}

// sign switchting correction
generated quantities{
  
  vector[P] lambdaMainC;
  vector[P] lambdaCrossC;
  corr_matrix[Q] PsiC; 
  
  PsiC = Psi;
  lambdaMainC = lambdaMain;
  lambdaCrossC = lambdaCross;

// factor 1 sign switching correction [p = 1, marker item]  
   if(lambdaMain[1] < 0){
     lambdaMainC[1:3] = -1*lambdaMain[1:3];
     lambdaCrossC[4:6] = -1*lambdaCross[4:6];

    if(lambdaMainC[4] > 0){ 
        PsiC[1, 2] = -1*Psi[1, 2];
        PsiC[2, 1] = -1*Psi[1, 2];
    }
}

// factor 2
   if(lambdaMain[4] < 0){
     lambdaMainC[4:6] = -1*lambdaMain[4:6];
    lambdaCrossC[1:3] = -1*lambdaCross[1:3];
     if(lambdaMain[1] > 0){
       PsiC[1, 2] = -1*Psi[1, 2];
       PsiC[2, 1] = -1*Psi[1, 2];
    }
  }
 
}
