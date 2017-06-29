
  real lognormal_mu(real mean, real variance) {
    return log(mean/(sqrt(1 + variance/(mean * mean))));
  }
  
  real lognormal_sd(real mean, real variance) {
    return sqrt(log(1 + variance/(mean * mean)));
  }