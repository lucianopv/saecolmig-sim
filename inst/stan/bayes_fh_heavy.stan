data {
  int<lower=1> N_in;
  int<lower=0> N_oos;
  int<lower=1> K;
  matrix[N_in, K]  X_in;
  matrix[N_oos, K] X_oos;
  vector[N_in]     y_direct;
  vector<lower=0>[N_in] vardir;
  real<lower=0>    sigma_v_scale;
}
parameters {
  vector[K]              beta;
  real<lower=0>          sigma_v;
  vector[N_in]           z;
}
transformed parameters {
  vector[N_in] v        = sigma_v * z;
  vector[N_in] theta_in = X_in * beta + v;
}
model {
  // Half-Normal on sigma_v (the <lower=0> constraint makes it half). Its
  // light tail, the covariates standardised in the R wrapper and
  // beta ~ N(0, 1) keep sigma_v identified.
  sigma_v ~ normal(0, sigma_v_scale);
  beta    ~ normal(0, 1);
  z       ~ normal(0, 1);
  y_direct ~ normal(theta_in, sqrt(vardir));
}
generated quantities {
  // OOS areas have no observed direct estimate, so their random effect
  // cannot be inferred. It is drawn from its prior N(0, sigma_v), so the
  // OOS posterior interval carries both the regression uncertainty and
  // the random-effect variance; X_oos * beta alone would omit the latter
  // and give intervals that are too narrow.
  vector[N_oos] theta_oos;
  for (i in 1:N_oos) {
    theta_oos[i] = dot_product(X_oos[i], beta) + normal_rng(0, sigma_v);
  }
}
