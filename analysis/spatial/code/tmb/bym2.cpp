#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() () {
  // Data
  DATA_VECTOR(y);       // Observed successes (binary or counts)
  DATA_VECTOR(N);       // Number of trials (for binomial)
  DATA_VECTOR(x);       // Covariate
  DATA_SPARSE_MATRIX(adj_matrix); // Adjacency matrix

  // Parameters
  PARAMETER(intercept); // Fixed intercept
  PARAMETER(HbAS_or_SS);      // Covariate coefficient
  PARAMETER_VECTOR(u);  // Unstructured random effects
  PARAMETER_VECTOR(v);  // Structured spatial effects (CAR)
  PARAMETER(log_tau_u); // Log-precision for unstructured effects
  PARAMETER(log_tau_v); // Log-precision for structured effects

  // In BYM formulation, parameter v must sum to 0 (to be identifiable)
  // Sp we pass in as n-1 dimensions and set nth value accordingly
  vector<Type> v2(v.size() + 1) ;
  v2.head(v.size()) = v ;  // Assign free elements
  v2(v.size()) = -v.sum();  // Enforce sum to zero

  // Transform precision parameters
  Type tau_u = exp(log_tau_u);
  Type tau_v = exp(log_tau_v);

  // Initialize negative log-likelihood
  Type nll = 0.0;

  // Priors for fixed effects (Bayesian sense)
  // Normal priors on intercept and beta
  nll -= dnorm(intercept, Type(0), Type(10), true);  // Mean 0, sd 10 for intercept (weakly informative)
  nll -= dnorm(HbAS_or_SS, Type(0), Type(10), true);       // Mean 0, sd 10 for beta (weakly informative)

  // Half-Cauchy prior on precisions (in log scale for positivity)
   nll -= dnorm(log_tau_u, Type(0), Type(1), true) - log(2);  // Half-Cauchy prior for tau_u
   nll -= dnorm(log_tau_v, Type(0), Type(1), true) - log(2);  // Half-Cauchy prior for tau_v

  // more informative priors option
  // nll -= dnorm(log_tau_u, Type(log(100)), Type(1), true) - log(2); // Expecting variance around 1/10 (precision of 10)
  // nll -= dnorm(log_tau_v, Type(log(100)), Type(1), true) - log(2);

  // Unstructured random effects (prior already included)
  nll -= sum(dnorm(u, Type(0), Type(1)/sqrt(tau_u), true));  // This is already a prior

  // Compute the Laplacian matrix Q = D - A
  // Construct Laplacian matrix (Q = D - A)//
  matrix<Type> Q =  - adj_matrix;  // Copy adjacency matrix
 
  // Compute row sums for sparse adjacency matrix
  Type penalty = 0.0 ;
  for (int k = 0; k < adj_matrix.outerSize(); ++k) {
    for (typename Eigen::SparseMatrix<Type>::InnerIterator it(adj_matrix, k); it; ++it) {
      // guard against counting the diagonal, in case
      // the adjacency matrix is encoded that way
      auto i = it.row() ;
      auto j = it.col() ;
      if( it.row() != it.col() && it.value() > 0 ) {
          penalty += -0.5 * (v2(i)-v2(j)) * (v2(i)-v2(j)) ;
          Q(it.row(), it.row()) += 1 ;
        }
    }
}

// Eg
// Q = [
//        -1  2  -1  0    0   0  0 0 0
//         0 -1   3  -1  -1  -1  0 0 0
// etc
// ]
// Then Qv =
// [ -v1 + 2v2 -v3      ]
// [ -v2 + 3v3 -v4 -v5  ]
//
// and v^t Q v
//
// [ -v1^2 + 2v1v2 -v1v3 -v2^2 +3v2v3 - v2v4 -v2v5 ... ]
//
// 

// ICAR penalty term
vector<Type> Qv = Q * v2;

// std::cerr << "nllbefore = " << nll << ".\n"; 
Type penalty2 = -Type(0.5) * (v2.transpose() * Qv).sum();

std::cerr << "v[1:4] = " << v[0] << ", " << v[1] << ", " << v[2] << ", " << v[3] << ".\n" ;
std::cerr << "penalty calculation: original = " << penalty2 << ", new = " << penalty << "\n" ;
nll -= penalty ;

// std::cerr << "nllafter = " << nll << ".\n"; 

//  std::cerr << "A[1,1:4] = " << adj_matrix.coeff(0,0) << " " << adj_matrix.coeff(0,1) << " " 
//            << adj_matrix.coeff(0,2) << " " << adj_matrix.coeff(0,3) << ".\n";

// // std::cerr << "D[1:4] = " << D(0) << " " << D(1) << " " 
// //           << D(2) << " " << D(3) << ".\n";

std::cerr << "Q[1,1:4] = " << Q(0,0) << " " << Q(0,1) << " " 
          << Q(0,2) << " " << Q(0,3) << ".\n";


// Not sure about this - BYM doesn't seem to have it
// Could be the scaling thing or the BYM2 thing
// Commenting out for now
// vector<Type> v_effect = Qv/tau_v;  // Doing the scaling as a post-step. 100% OK!!

  // Linear predictor and binomial likelihood
  vector<Type> eta = intercept + HbAS_or_SS * x + u + v2 ; // not sure about QV/tau here._effect;
  vector<Type> p = invlogit(eta); // Logit link for binomial
  nll -= sum(dbinom(y, N, p, true)); 

// std::cerr << "nllfinal = " << nll << ".\n"; 

  return nll;
}
