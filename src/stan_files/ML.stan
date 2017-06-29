
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
  
  vector<lower=0>[nbreaks+1] muZ;
  vector<lower=0>[nbreaks+1] sdZ;
  vector<lower=0>[nbreaks+1] alpha_dirichlet;
  real<lower=0> musigma;
  real<lower=0> sdsigma; 
}

parameters {
  vector<lower=0>[nbreaks+1] Z;
  simplex[nbreaks+1] Z_duration;
  real<lower=0> sigma;
}

transformed parameters {
  vector<lower=0>[nbreaks] yearZ;
  yearZ = cumulative_sum(Z_duration[1:nbreaks]) * count;
}

model {
  vector[count] Lpred;
  
  // Priors for Z, yearZ, sigma
  for (i in 1:(nbreaks+1)) Z[i] ~ lognormal(lognormal_mu(muZ[i], square(sdZ[i])), lognormal_sd(muZ[i], square(sdZ[i]))) T[0, 3];
  Z_duration ~ dirichlet(alpha_dirichlet);
  sigma ~ lognormal(lognormal_mu(musigma, square(sdsigma)), lognormal_sd(musigma, square(sdsigma)));
  
  Lpred = generate_Lpred(nbreaks, count, Lobs, ss, Lc, Linf, K, Z, yearZ);
  for(m in 1:count) {
    if(ss[m]>0) target += normal_lpdf(Lobs[m] | Lpred[m], sigma/sqrt(ss[m]));
	//if(ss[m]>0) Lobs[m] ~ normal(Lpred[m], sigma/sqrt(ss[m]));
  }
}

generated quantities {
  vector[count] Lpred;
  vector[count] log_lik;
  Lpred = generate_Lpred(nbreaks, count, Lobs, ss, Lc, Linf, K, Z, yearZ);
  for (m in 1:count) log_lik[m] = ss[m] > 0 ? normal_lpdf(Lobs[m] | Lpred[m], sigma/sqrt(ss[m])) : 0;
}

