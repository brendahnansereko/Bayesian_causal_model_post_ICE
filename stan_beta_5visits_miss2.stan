data {
  int<lower=1> K; // Number of time points
  int<lower=1> J; // Number of covariates
  int<lower=0> N; // Number of units
  vector[J] x[N]; // Covariates
  vector[K] y[N]; // Outcome matrix
  int<lower=0, upper=1> W[N]; // Treatment indicator
  int<lower=0> num_discont; // Number of discontinuators
  int discont_index[num_discont]; // Indices of discontinuators
  int Z[N]; // Group indicator
  int numMissPer0[K]; // Missing counts for group 0
  int numMissPer1[K]; // Missing counts for group 1
  vector[K] priorMiss; // Prior for missingness
  int<lower=1, upper=K> t_indiv[N]; // Individual discontinuation times
  int<lower=0> num_discontonly; // Number of discontinuities only no missing
  int discontonly_index[num_discontonly]; // discontinuities only no missing
}

parameters {
  matrix[K, J] beta; // Covariate effects
  cholesky_factor_corr[K] L; // Cholesky factor of correlation
  vector<lower=0>[K] tau; // Scale
  vector[K] alpha0; // Control group mean
  vector[K] alpha1; // Treated group mean
  real k_0; // Causal parameter
  simplex[K] propMisPer0;
  simplex[K] propMisPer1;
}

transformed parameters {
  cov_matrix[K] Sigma;
  Sigma = quad_form_diag(L * L', tau); // Full covariance
 // real k_0 = 1;

}

model {
  // Priors
  to_vector(beta) ~ normal(0, 100);
  alpha0 ~ normal(0, 100);
  alpha1 ~ normal(0, 100);
  k_0 ~ normal(1, 1);
  propMisPer0 ~ dirichlet(priorMiss);
  propMisPer1 ~ dirichlet(priorMiss);

  // Model for non-discontinuators
  for (n in 1:N) {
  int t_n = t_indiv[n] ;
    vector[t_n] mu_pre = W[n] * alpha1[1:t_n] + (1 - W[n]) * alpha0[1:t_n] + beta[1:t_n] * x[n];
    matrix[t_n, t_n] Sigma_pre = Sigma[1:t_n, 1:t_n];
    y[n][1:t_n] ~ multi_normal(mu_pre, Sigma_pre);
  }

  

// Model for discontinuators
for (z in 1:num_discontonly) {
  int n = discontonly_index[z];      // Index of discontinuing subject
  int t_n = t_indiv[n];       // Time point of discontinuation (in Stan, time indices often start from 1)
  int post_n = K - t_n;           // Number of post-discontinuation points

    // Covariance submatrices
    matrix[t_n, t_n] Sigma_11 = Sigma[1:t_n, 1:t_n];
    matrix[post_n, t_n] Sigma_21 = Sigma[(t_n + 1):K, 1:t_n];
    matrix[t_n, post_n] Sigma_12 = transpose(Sigma_21);
    matrix[post_n, post_n] Sigma_22 = Sigma[(t_n + 1):K, (t_n + 1):K];
    matrix[t_n, t_n] Sigma_11_inv = inverse(Sigma_11);
    matrix[post_n, post_n] Sigma_post = Sigma_22 - Sigma_21 * Sigma_11_inv * Sigma_12;
    matrix[post_n, t_n] beta_t = Sigma_21 * Sigma_11_inv;

    // Conditional means
    vector[post_n] mu_post = (beta_t * (y[n][1:t_n] - (alpha1[1:t_n] + beta[1:t_n]*x[n]))) 
                        + rep_vector(k_0 * (alpha1[t_n] - alpha0[t_n]), post_n) 
                        + (alpha0[(t_n + 1):K] + beta[(t_n + 1):K] * x[n]);

    y[n][(t_n + 1):K] ~ multi_normal(mu_post, Sigma_post);
  }


  // Missingness model
  numMissPer0 ~ multinomial(propMisPer0);
  numMissPer1 ~ multinomial(propMisPer1);
}

generated quantities {
  real MAR = alpha1[K] - alpha0[K];
  real J2R = propMisPer1[K] * (alpha1[K] - alpha0[K]);
  real CIR = sum(propMisPer1[1:K] .* (alpha1[1:K] - alpha0[1:K]));
  real causal_J2R = J2R + k_0 * sum(propMisPer1[1:(K-1)] .* (alpha1[1:(K-1)] - alpha0[1:(K-1)]));
}
