#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() () {
  // Data
  DATA_VECTOR(y);       // Observed successes (binary or counts)
  DATA_VECTOR(N);       // Number of trials (for binomial)
  DATA_VECTOR(x);       // Covariate
  
  // Parameters
  PARAMETER(intercept); // Fixed intercept
  PARAMETER(HbAS_or_SS);      // Covariate coefficient

  // Initialize negative log-likelihood
  Type nll = 0.0;

  // Priors for fixed effects (Bayesian sense)
  // Normal priors on intercept and HbAS_or_SS
  nll -= dnorm(intercept, Type(0), Type(10), true);  // Mean 0, sd 10 for intercept (weakly informative)
  nll -= dnorm(HbAS_or_SS, Type(0), Type(10), true);       // Mean 0, sd 10 for beta of HbAS_or_SS (weakly informative)

  // Linear predictor and binomial likelihood
  vector<Type> eta = intercept + HbAS_or_SS * x;
  vector<Type> p = invlogit(eta); // Logit link for binomial
  nll -= sum(dbinom(y, N, p, true)); 

  return nll;
}
