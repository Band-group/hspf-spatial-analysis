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
}

template< typename Type >
Type objective_function<Type>::operator() () {
	// --- Data ---
	DATA_VECTOR(y);       // Observed successes (binary or counts)
	DATA_VECTOR(N);       // Number of trials (for binomial)
	DATA_VECTOR(x);       // Covariate

	// --- Parameters ---
	PARAMETER(        intercept   ) ; // Fixed intercept
	PARAMETER(        beta        ) ; // Covariate coefficient
	PARAMETER_VECTOR( u           ) ; // Unstructured random effects
	PARAMETER_VECTOR( v           ) ; // Structured spatial effects (CAR)
	PARAMETER(        log_tau_u   ) ; // Log-precision for unstructured effects
	PARAMETER(        log_tau_v   ) ; // Log-precision for structured effects
	// Q matrix
	// This should be derived from the adjacency matrix A as
	// Q_ii = number of neighbours of cell i (not including i)
	// Q_ij = -1 if i and j are neighbours (i!=j)
	//
	// Optionally Q can also be 'rescaled' 
	DATA_SPARSE_MATRIX(Q); 

	// --- Transforms ---
	// Structured effects:
	// In BYM formulation, parameter v must sum to 0 (to be identifiable)
	// So we pass in as n-1 dimensions and set nth value accordingly
	assert( v.size() == Q.rows() - 1 ) ;
	vector<Type> va = append_value( v , -v.sum() ) ;

	// Transformed precision parameters
	Type tau_u = exp(log_tau_u);
	Type tau_v = exp(log_tau_v);

	// --- LL calculation ---
	Type nll = 0.0; // negative log-likelihood

	nll -= log_dN( intercept,  0.0, 10.0 ) ;		// Mean 0, sd 10 for intercept (weakly informative)
	nll -= log_dN( beta, 0.0, 10.0 ) ;				// Mean 0, sd 10 for beta (weakly informative)

	nll -= log_dN( log_tau_u, 0.0, 1.0 ) + log(2) ;		// half-normal on log tau_u
	nll -= log_dN( log_tau_v, 0.0, 1.0 ) + log(2) ;		// half-normal on log tau_v

	nll -= sum( log_dN( u, Type(0.0), 1.0/sqrt(tau_u) )) ; // independent random effects are gaussian with specified precision

	// ICAR (spatial) penalty term
	vector<Type> Qv = Q * va;
	Type penalty = -Type(0.5) * (va * Qv).sum() ;
	nll -= penalty ;

	std::cerr << "Q[1,1:4] = " << Q.coeff(0,0) << " " << Q.coeff(0,1) << " " << Q.coeff(0,2) << " " << Q.coeff(0,3) << ".\n";
	std::cerr << "v[1:4] = " << v[0] << ", " << v[1] << ", " << v[2] << ", " << v[3] << ".\n" ;
	std::cerr << "Qv[1:4,156,158] = " << Qv[0] << ", " << Qv[1] << ", " << Qv[2] << ", " << Qv[3] << ", " << Qv[155] << ", " << Qv[157] << ".\n" ;
	std::cerr << "penalty: " << penalty << "\n" ;

	// Not sure about this - BYM doesn't seem to have it
	// Could be the scaling thing or the BYM2 thing
	// Commenting out for now
	// vector<Type> v_effect = Qv/tau_v;  // Doing the scaling as a post-step. 100% OK!!

	// Linear predictor and binomial likelihood
	vector<Type> eta = intercept + beta*x + u + va ; // not sure about QV/tau here._effect;
	vector<Type> p = invlogit(eta); // Logit link for binomial
	nll -= sum(dbinom(y, N, p, true)); 

	// std::cerr << "nllfinal = " << nll << ".\n"; 

	return nll;
}
