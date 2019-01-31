
functions {
#include /functions/misc.stan
#include /functions/generate_Lpred.stan
}

data {
  int<lower=0> count;
  int<lower=1> nbreaks;
  vector<lower=0>[count] Lobs;
  vector<lower=0>[count] ss;
  real<lower=0> Lc;  
  real<lower=0> Linf;
  real<lower=0> K;
  
  int<lower=0, upper=1> Z_dist;
  matrix<lower=0>[nbreaks+1, 2] Z_par;
  vector<lower=0>[nbreaks+1] alpha_dirichlet;
  int<lower=0, upper=1> sigma_dist;
  vector<lower=0>[2] sigma_par;
}

transformed data {
  matrix<lower=0>[nbreaks+1, 2] Z_par_ln;
  vector<lower=0>[2] sigma_par_ln;
  
  for (i in 1:(nbreaks+1)) {
    Z_par_ln[i, 1] = Z_dist == 1 ? lognormal_mu(Z_par[i, 1], Z_par[i, 2]) : 0;
    Z_par_ln[i, 2] = Z_dist == 1 ? lognormal_sd(Z_par[i, 1], Z_par[i, 2]) : 0;
  }
  
  sigma_par_ln[1] = sigma_dist == 1 ? lognormal_mu(sigma_par[1], sigma_par[2]) : 0;
  sigma_par_ln[2] = sigma_dist == 1 ? lognormal_sd(sigma_par[1], sigma_par[2]) : 0;
}

parameters {
  vector<lower=0>[nbreaks+1] Z;
  simplex[nbreaks+1] Z_duration;
  real<lower=0> sigma;
}

transformed parameters {
  vector<lower=0>[nbreaks] D;
  D = cumulative_sum(Z_duration[1:nbreaks]) * (count - 1) + 1;
}

model {
  vector[count] Lpred;
  
  // Priors
  for (i in 1:(nbreaks+1)) {
    if (Z_dist == 0) Z[i] ~ uniform(Z_par[i, 1], Z_par[i, 2]);
    if (Z_dist == 1) Z[i] ~ lognormal(Z_par_ln[i, 1], Z_par_ln[i, 2]) T[0, 3];
  }
  
  Z_duration ~ dirichlet(alpha_dirichlet);
  
  if (sigma_dist == 0) sigma ~ uniform(sigma_par[1], sigma_par[2]);
  if (sigma_dist == 1) sigma ~ lognormal(sigma_par_ln[1], sigma_par_ln[2]);
  
  // Generate predicted mean length 
  Lpred = generate_Lpred(nbreaks, count, Lobs, ss, Lc, Linf, K, Z, D);
  
  // Likelihood
  for (m in 1:count) {
    if (ss[m]>0) target += normal_lpdf(Lobs[m] | Lpred[m], sigma/sqrt(ss[m]));
  }
}

generated quantities {
  vector[count] Lpred;
  vector[count] log_lik;
  Lpred = generate_Lpred(nbreaks, count, Lobs, ss, Lc, Linf, K, Z, D);
  for (m in 1:count) log_lik[m] = ss[m] > 0 ? normal_lpdf(Lobs[m] | Lpred[m], sigma/sqrt(ss[m])) : 0;
}

