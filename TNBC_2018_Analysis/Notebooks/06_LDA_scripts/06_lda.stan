// copyright to Kris Sankaran https://github.com/krisrs1128/boot_expers/blob/master/ldaSim/inst/extdata/lda.stan.
data {
  int<lower=1> K; // num topics
  int<lower=1> V; // num words (ASVs)
  int<lower=1> D; // num docs (samples)
  int<lower=0> n[D, V]; // word (ASVs) counts for each doc (sample)
  // n is the x; D is rows, V is the columns

  // hyperparameters
  vector<lower=0>[K] alpha;
  vector<lower=0>[V] gamma;
}

parameters {
  simplex[K] theta[D]; // topic mixtures
  simplex[V] beta[K]; // word (ASVs) dist for k^th topic
}

model {
  for (d in 1:D) {
    theta[d] ~ dirichlet(alpha);
  }

  for (k in 1:K) {
    beta[k] ~ dirichlet(gamma);
  }

  for (d in 1:D) {
    vector[V] eta;
    eta = beta[1] * theta[d, 1];

    for (k in 2:K) {
      eta = eta + beta[k] * theta[d, k];
    }
    
    n[d] ~ multinomial(eta);
    
  }
}

generated quantities {
  int<lower=0> x_sim[D, V]; // simulated word counts, for posterior checking
  vector[D] log_lik;// log conditional likelihood
  
  for (d in 1:D) {
    vector[V] eta;
    eta = beta[1] * theta[d, 1];

    for (k in 2:K) {
      eta = eta + beta[k] * theta[d, k];
    }
    
    x_sim[d] = multinomial_rng(eta, sum(n[d]));
    log_lik[d] = multinomial_lpmf(n[d] | eta);
  }
}
