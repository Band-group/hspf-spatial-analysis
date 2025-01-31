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
  PARAMETER_VECTOR(u);  // Unstructured random effects
  PARAMETER(log_tau_u); // Log-precision for unstructured effects
  
  // Transform precision parameters
  Type tau_u = exp(log_tau_u);
  
  // Initialize negative log-likelihood
  Type nll = 0.0;

  // Priors for fixed effects (Bayesian sense)
  // Normal priors on intercept and HbAS_or_SS
  nll -= dnorm(intercept, Type(0), Type(10), true);  // Mean 0, sd 10 for intercept (weakly informative)
  nll -= dnorm(HbAS_or_SS, Type(0), Type(10), true);       // Mean 0, sd 10 for beta of HbAS_or_SS (weakly informative)

  // Half-Cauchy prior on precisions (in log scale for positivity)
  // nll -= dnorm(log_tau_u, Type(0), Type(1), true) - log(2);  // Half-Cauchy prior for tau_u
  
  // more informative priors option
  nll -= dnorm(log_tau_u, Type(log(100)), Type(1), true) - log(2); // Expecting variance around 1/10 (precision of 10)
  
  // Unstructured random effects (prior already included)
  nll -= sum(dnorm(u, Type(0), Type(1)/sqrt(tau_u), true));  // This is already a prior

  // Linear predictor and binomial likelihood
  vector<Type> eta = intercept + HbAS_or_SS * x + u;
  vector<Type> p = invlogit(eta); // Logit link for binomial
  nll -= sum(dbinom(y, N, p, true)); 

  return nll;
}
