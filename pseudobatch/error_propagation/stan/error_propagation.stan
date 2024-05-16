/* A Stan model for error propagation in pseudo batch transformation.

   This model assumes that the errors for measurements of compartment volume,
   species concentrations, sample volumes and feed concentrations are all known.

 */
functions {
#include functions.stan
}
data {
  int<lower=2> N;   // number of observations
  int<lower=1> S;   // number of species
  vector[N] y_v;
  array[N, S] real y_c;
  vector[N] y_s;
  vector[N] y_f;
  vector[S] y_cfeed;
  real<lower=0> sigma_v;
  vector<lower=0>[S] sigma_c;
  real<lower=0> sigma_s;
  real<lower=0> sigma_f;
  real<lower=0> sigma_cfeed;
  array[2] real prior_apump;
  array[2] real prior_as;
  array[2] real prior_v0;
  array[2, S] real prior_m;
  array[2] real prior_f_nonzero;
  array[2, S] real prior_cfeed;
  int<lower=0, upper=1> likelihood;
}
transformed data {
  int N_f_zero = count_zeros(y_f);
  int N_s_zero = count_zeros(y_s);
  int N_cfeed_nonzero = S - count_zeros(y_cfeed);
  int N_f_nonzero = N - N_f_zero;
  int N_s_nonzero = N - N_s_zero;
  array[N_f_zero] int ix_f_zero = vector_zero_ixs(y_f);
  array[N_f_nonzero] int ix_f_nonzero = vector_nonzero_ixs(y_f);
  array[N_s_nonzero] int ix_s_nonzero = vector_nonzero_ixs(y_s);
  array[N_cfeed_nonzero] int ix_cfeed_nonzero = vector_nonzero_ixs(y_cfeed);
}
parameters {
  real<lower=0> v0; // starting volume
  matrix<lower=0>[N, S] m;  // mass of each species
  vector[N] as; // logit fraction of current volume sampled at each point
  vector<lower=0>[N_f_nonzero] f_nonzero; // amount of feed in interval prior to each point
  vector<lower=0>[N_cfeed_nonzero] cfeed_nonzero; // concentration of feed
  real apump; // log-scale factor by which the pump is biased
}
transformed parameters {
  vector[N] f;
  vector[N] s;
  vector[N] v;
  matrix[N, S] c;
  f[ix_f_zero] = zeros_vector(N_f_zero);
  f[ix_f_nonzero] = f_nonzero;
  v[1] = v0 + f[1];
  s[1] = v[1] * inv_logit(as[1]);
  for (n in 2 : N) {
    v[n] = v0 + sum(f[ : n]) - sum(s[ : n - 1]);
    s[n] = v[n] * inv_logit(as[n]);
  }
  for (species_i in 1 : S){
    c[,species_i] = m[,species_i] ./ v;
  }
}
model {
  apump ~ normal(prior_apump[1], prior_apump[2]);
  as ~ normal(prior_as[1], prior_as[2]);
  v0 ~ lognormal(prior_v0[1], prior_v0[2]);
  for (species_i in 1 : S){
    m[,species_i] ~ lognormal(prior_m[1, species_i], prior_m[2, species_i]);
  }
  f_nonzero ~ lognormal(prior_f_nonzero[1], prior_f_nonzero[2]);
  if (likelihood) {
    y_v ~ lognormal(log(v), sigma_v);
    for (species_i in 1 : S){
      y_c[,species_i] ~ lognormal(log(c[,species_i]), sigma_c[species_i]);
    }
    y_s[ix_s_nonzero] ~ lognormal(log(s[ix_s_nonzero]), sigma_s);
    y_f[ix_f_nonzero] ~ lognormal(log(f[ix_f_nonzero]) + apump, sigma_f);
  }
  if (N_cfeed_nonzero > 0){
    cfeed_nonzero ~ lognormal(prior_cfeed[1, ix_cfeed_nonzero],
                              prior_cfeed[2, ix_cfeed_nonzero]);
    if (likelihood){
      y_cfeed[ix_cfeed_nonzero] ~ lognormal(log(cfeed_nonzero), sigma_cfeed);
    }
  }
}
generated quantities {
  real pump_bias = log(apump);
  vector[S] cfeed = zeros_vector(S);
  matrix[N, S] pseudobatch_c;
  cfeed[ix_cfeed_nonzero] = cfeed_nonzero;
  for (species_i in 1 : S){
    pseudobatch_c[, species_i] = pseudobatch_transform(v, 
                                                       s, 
                                                       c[,species_i], 
                                                       f, 
                                                       cfeed[species_i]);
  }
}

