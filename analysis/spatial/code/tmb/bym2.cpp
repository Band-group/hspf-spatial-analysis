#include <TMB.hpp>
#include <cassert>

// #define DEBUG 1

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
	DATA_SCALAR( prior_intercept_sd );
	DATA_SCALAR( prior_beta_sd );
	DATA_SCALAR( prior_sd_rate );

	// --- Parameters ---
	PARAMETER(         intercept ) ; // Fixed intercept
	PARAMETER(              beta ) ; // Covariate coefficient
	PARAMETER_VECTOR(          u ) ; // Unstructured random effects
	PARAMETER_VECTOR(          v ) ; // Structured spatial effects (CAR)
	PARAMETER(           log_tau ) ; // Log-precision for random effects. Passed as log(tau) so that tau is enforced positive.
	PARAMETER(       logodds_phi ) ; // Log-odds of mixture proportion of spatial vs independent components.
	// Q matrix
	// This should be derived from the adjacency matrix A as
	// Q_ii = number of neighbours of cell i (not including i)
	// Q_ij = -1 if i and j are neighbours (i!=j)
	//
	// Optionally Q can also be 'rescaled' 
	DATA_SPARSE_MATRIX(Q); 
	DATA_MATRIX(connected_components); 

	// Transformed precision parameters
	Type tau = exp(log_tau) ;
	Type sd_of_random_effects = 1/sqrt(tau) ;
	Type phi = exp(logodds_phi) / ( 1 + exp(logodds_phi)) ;

#if DEBUG
	std::cerr << "params:\n"
		<< " -     log(tau): " << log_tau << "\n"
		<< " -          tau: " << tau << "\n"
		<< " - logodds(phi): " << logodds_phi << "\n"
		<< " -          phi: " << phi << "\n"
		<< " -         beta: " << beta << "\n"
		<< " -    intercept: " << intercept << ".\n" ;
	std::cerr
		<< " -     Q[1,1:4]: " << Q.coeff(0,0) << " " << Q.coeff(0,1) << " " << Q.coeff(0,2) << " " << Q.coeff(0,3) << ".\n";
#endif

	// --- Transforms ---

	// --- LL calculation ---
	Type nll = 0.0; // negative log-likelihood

	// Spatial coefficient v is supposed to sum to 0
	// Solve here by requiring a very narrow normal prior on the sum
	// of values in each connected component.
	nll -= sum(
		log_dN(
			connected_components * v,
			0.0,
			0.00001
		)
	) ;

	// Priors
	nll -= log_dN( intercept,  0.0, prior_intercept_sd ) ;		 // Weak 0-centred prior on intercept
	nll -= log_dN( beta,       0.0, prior_beta_sd      ) ;		 // Weak 0-centred prior on beta

	// Prior on random effect sd
	// Riebler et al (2016) says this is a "Type 2 Gumbel" distribution on tau,
	// or a an exponential prior on the sd:
	nll -= dexp( sd_of_random_effects, prior_sd_rate, /* give_log */ 1 ) ;
	nll -= sum( log_dN( u, 0.0, 1.0 )) ; // independent random effects are standard gaussian

	// ICAR (spatial) penalty term
	vector<Type> Qv = Q * v;
	Type penalty = -Type(0.5) * (v * Qv).sum() ;
	nll -= penalty ;

	// Linear predictor and binomial likelihood
	vector<Type> predictor = (
		intercept
		+ (beta * x)
		+ sd_of_random_effects * sqrt(1-phi)*u + sd_of_random_effects * sqrt(phi)*v
	) ;
	vector<Type> p = invlogit( predictor ) ; // Logit link for binomial
	nll -= sum( dbinom(y, N, p, true) ); 

	return nll;
}
