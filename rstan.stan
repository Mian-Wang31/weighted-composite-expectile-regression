data {
  int<lower=0> N;
  int<lower=0> p;
  vector[N] y;
  matrix[N,p] x;
  int<lower=0> L;
  vector[L] tau;
}

parameters {
  vector[L] b0;
  vector[p] beta;
  real<lower=0> sigma2;
  vector<lower=0>[p] delta2;
  simplex[L] w;
}


model {
  for (i in 1:p) {
  beta[i] ~ normal(0,delta2[i]);
  delta2[i] ~ inv_gamma(0.5,0.5);
  }
  sigma2 ~ inv_gamma(0.0001,1);
  w ~ dirichlet(rep_vector(0.1,L));
  for(i in 1:N){
    for(l in 1:L){
      real mu=b0[l]+x[i]*beta;
      if(y[i]<=mu){
        target += log(w[l]*((2/sqrt(sigma2*pi()))*((sqrt(1/(1-tau[l]))+sqrt(1/tau[l]))^(-1))*exp(-((1-tau[l])*(y[i]-mu)^2)/sigma2)));
      } else {
         target += log(w[l]*((2/sqrt(sigma2*pi()))*((sqrt(1/(1-tau[l]))+sqrt(1/tau[l]))^(-1))*exp(-(tau[l]*(y[i]-mu)^2)/sigma2)));
      }
    }
  }
  
}
