#include <TMB.hpp>
#include <cassert>

namespace {
	// Helper functions
	// These go in an un-named namespace, because in C++ this ensures they don't conflict with
	// things defined elsewhere.
	template< typename Type >
	vector< Type > append_value( vector<Type> v, Type value ) {
		vector< Type > result( v.size() + 1 ) ;
		result.head( v.size() ) = v ;
		result[ v.size() ] = value ;
		return result ;
	}

	// Helper to make writing log normal densities slicker.
	// Parameter types are a bit complex because it is used with 
	// doubles, Types, and vectors of Types as well.
	template< typename T1, typename T2, typename T3 >
	T1 log_dN( T1 value, T2 mean, T3 sd ) {
		return dnorm( value, T1(mean), T1(sd), 1 ) ;
	}

	// Version of the above with vector first parameter and returning a vector
	template< typename T1, typename T2, typename T3 >
	tmbutils::vector<T1> log_dN( tmbutils::vector<T1> value, T2 mean, T3 sd ) {
		return dnorm( value, T1(mean), T1(sd), 1 ) ;
	}

	template< typename T1, typename T2, typename T3 >
	T1 log_dhalfN( T1 value, T2 mean, T3 sd ) {
		return(
			(value >= mean)
			?
			( log(2) + log_dN( value, mean, sd ) )
			: 
			-10000
		) ;
	}
}

template< typename Type >
Type objective_function<Type>::operator() () {
	// --- Data ---
	DATA_VECTOR(y);       // Observed successes (binary or counts)
	DATA_VECTOR(N);       // Number of trials (for binomial)
	DATA_VECTOR(x);       // Covariate

	// -- Prior parameters ---
	DATA_SCALAR( prior_halfnormal_sd_tau );
	DATA_SCALAR( prior_independent_sd );
	DATA_SCALAR( prior_intercept_sd );
	DATA_SCALAR( prior_beta_sd );

	// --- Parameters ---
	PARAMETER(         intercept ) ; // Fixed intercept
	PARAMETER(              beta ) ; // Covariate coefficient
	PARAMETER_VECTOR(          u ) ; // Unstructured random effects
	PARAMETER_VECTOR(          v ) ; // Structured spatial effects (CAR)
	PARAMETER(           log_tau ) ; // Log-precision for random effects
	PARAMETER(       logodds_phi ) ; // Log-odds of mixture proportion of spatial vs independent components.
	// Q matrix
	// This should be derived from the adjacency matrix A as
	// Q_ii = number of neighbours of cell i (not including i)
	// Q_ij = -1 if i and j are neighbours (i!=j)
	//
	// Optionally Q can also be 'rescaled' 
	DATA_SPARSE_MATRIX(Q); 

	// Transformed precision parameters
	Type tau = exp(log_tau) ;
	Type phi = exp(logodds_phi) / ( 1 + exp(logodds_phi)) ;

	std::cerr << "params:\n"
		<< " -     log(tau): " << log_tau << "\n"
		<< " -          tau: " << tau << "\n"
		<< " - logodds(phi): " << logodds_phi << "\n"
		<< " -          phi: " << phi << "\n"
		<< " -         beta: " << beta << "\n"
		<< " -    intercept: " << intercept << ".\n" ;
	std::cerr
		<< " -     Q[1,1:4]: " << Q.coeff(0,0) << " " << Q.coeff(0,1) << " " << Q.coeff(0,2) << " " << Q.coeff(0,3) << ".\n";

	// --- Transforms ---

	// --- LL calculation ---
	Type nll = 0.0; // negative log-likelihood

	// Spatial coefficient v is supposed to sum to 0
	// (Actually on each connected component - multiple ccs not implemented yet)
	// Implement this here by having a strict prior on mean(v)
	nll -= log_dN( sum(v)/v.size(), 0.0, 0.00001 ) ;

	// Priors
	nll -= log_dN( intercept,  0.0, prior_intercept_sd ) ;		 // Weak 0-centred prior on intercept
	nll -= log_dN( beta,       0.0, prior_beta_sd      ) ;		 // Weak 0-centred prior on beta
	// half-normal prior on tau
	nll -= log_dhalfN( tau, 0.0, prior_halfnormal_sd_tau ) ; // normal on log of precision
	nll -= sum( log_dN( u, 0.0, prior_independent_sd )) ; // independent random effects are standard gaussian

	// ICAR (spatial) penalty term
	vector<Type> Qv = Q * v;
	Type penalty = -Type(0.5) * (v * Qv).sum() ;
	nll -= penalty ;

	std::cerr << "v[1:4] = " << v[0] << ", " << v[1] << ", " << v[2] << ", " << v[3] << ".\n" ;
	std::cerr << "Qv[1:4,156,158] = " << Qv[0] << ", " << Qv[1] << ", " << Qv[2] << ", " << Qv[3] << ", " << Qv[155] << ", " << Qv[157] << ".\n" ;
	std::cerr << "penalty: " << penalty << "\n" ;

	// Linear predictor and binomial likelihood
	Type sd_of_random_effects = 1/sqrt(tau) ;
	vector<Type> predictor = (
		intercept
		+ (beta * x)
		//+ sd_of_random_effects * sqrt(1-phi)*u + sd_of_random_effects * sqrt(phi)*v
	) ;
	vector<Type> p = invlogit( predictor ) ; // Logit link for binomial
	nll -= sum( dbinom(y, N, p, true) ); 

	return nll;
}
