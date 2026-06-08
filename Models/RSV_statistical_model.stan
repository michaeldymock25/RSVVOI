// MCMC implementation of the RSV statistical model
// strategies are NS (d == 1), MV (d == 2) and mAbs (d == 3)

data {
  array[3] int n;                 // number of participants randomised to each strategy
  array[3] int y;                 // number of events for each strategy
  real<lower=0,upper=1> alpha;    // power prior contribution
  real a_NS;                      // first shape hyperparameter for NS event probability
  real b_NS;                      // second shape hyperparameter for NS event probability
  vector[2] mu;                   // mean hyperparameters for log relative effect parameters
  vector<lower=0>[2] sigma;       // standard deviation hyperparameters for log relative effect parameters
  int<lower=0,upper=1> OR;        // OR == 1 for odds ratio model or OR == 0 for relative risk model
}

parameters {
  real<lower=0,upper=1> p_NS;     // event probability for the NS strategy
  vector[2] pi;                   // log relative effects 
} 

transformed parameters {
  vector<lower=0,upper=1>[2] p_d; // event probability for the MV and mAbs strategies
  if(OR == 0) p_d = p_NS*exp(pi);
  if(OR == 1) p_d = inv_logit(logit(p_NS) + pi); 
}

model {
  vector[3] p = append_row([p_NS]', p_d);
  vector[2] sd;
  if(alpha == 0){
    sd = to_vector({1,1});
  } else {
    sd = sqrt(sigma^2/sqrt(alpha));
  }
  target += beta_lpdf(p_NS | 1 + alpha*(a_NS - 1), 1 + alpha*(b_NS - 1));
  target += normal_lpdf(pi | mu, sd);
  target += binomial_lpmf(y | n, p);
}

generated quantities {
  vector[2] r = p_d/p_NS;
}
