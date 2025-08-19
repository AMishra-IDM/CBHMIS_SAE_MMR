
functions {
  real standard_icar_lpdf(vector phi, array[,] int adj, real epsilon) {
    if (size(adj) != 2)
      reject("require 2 rows for adjacency array;",
             " found rows = ", size(adj));
    return -0.5 * dot_self(phi[adj[1]] - phi[adj[2]])
           + normal_lupdf(sum(phi) | 0, epsilon * rows(phi));
  }
}
data {
  int<lower=1> N;                      // number of regions
  int<lower=0> z[N];                   // observed maternal deaths
  vector<lower=0>[N] live_births;      // offset
  int<lower=1> K_x;                    // incidence covariates
  matrix[N, K_x] x;
  int<lower=1> K_w;                    // reporting covariates
  matrix[N, K_w] w;

  int<lower=0> N_edges;
  array[2, N_edges] int<lower=1, upper=N> neighbors;

  real<lower=0> tau; // scaling factor for ICAR
}
transformed data {
  vector[N] log_E = log(live_births);
}
parameters {
  real alpha0;
  vector[K_x] alpha;
  real beta0;
  vector[K_w] beta;

  real<lower=0, upper=1> rho;
  real<lower=0> sigma;

  vector[N] phi;
  vector[N] theta;
}
transformed parameters {
  vector[N] gamma;
  gamma = sigma * (sqrt(1 - rho) * theta + sqrt(rho / tau) * phi);
}
model {
  vector[N] log_lambda;
  vector[N] logit_pi;
  vector[N] log_mu;

  // Priors
  alpha0 ~ normal(0, 5);
  alpha ~ normal(0, 2);
  beta0 ~ normal(2, 0.6);  // informative prior on reporting rate
  beta ~ normal(0, 2);
  sigma ~ normal(0, 2);
  rho ~ beta(0.5, 0.5);

  phi ~ standard_icar(neighbors, 0.001);
  theta ~ std_normal();

  // Model
  log_lambda = alpha0 + x * alpha + gamma;
  logit_pi = beta0 + w * beta;
  log_mu = log_E + log_lambda + log_inv_logit(logit_pi);

  z ~ poisson_log(log_mu);
}
generated quantities {
  vector[N] lambda = exp(alpha0 + x * alpha + gamma);
  vector[N] pi = inv_logit(beta0 + w * beta);
  vector[N] mu = live_births .* lambda .* pi;

  array[N] int z_rep;
  for (i in 1:N) {
    z_rep[i] = poisson_log_rng(log(mu[i]));
  }
}
