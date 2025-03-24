#include <TMB.hpp>
#include <cassert>

#define DEBUG 1

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

	template< typename T1, typename T2, typename T3 >
	tmbutils::vector<T1> clamp( tmbutils::vector<T1> v, T2 min, T3 max ) {
		tmbutils::vector<T1> result = v ;
		for( int i = 0; i < v.size(); ++i ) {
			result(i) = std::max( std::min( result(i), T1(max) ), T1(min) ) ;
		}
		return result ;
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
	DATA_SCALAR( prior_logodds_phi_mean );
	DATA_SCALAR( prior_logodds_phi_sd );
	DATA_SCALAR( prior_log_nu_sd ) ;

	DATA_STRING( model_choice ); // "norandom" or "bym2"
	DATA_STRING( link_choice ); // "logit" or "generalised-logit" or "linear"
	Type NaN = 0.0/0.0 ;

	if(
		model_choice != "bym2"
		&& model_choice != "norandom"
	) {
		return NaN ;
	}

	if(
		link_choice != "logit"
		&& link_choice != "generalised-logit"
		&& link_choice != "linear"
	) {
		return NaN ;
	}

	// --- Parameters ---
	PARAMETER(         intercept ) ; // Fixed intercept
	PARAMETER(              beta ) ; // Covariate coefficient
	PARAMETER_VECTOR(          u ) ; // Unstructured random effects
	PARAMETER_VECTOR(          v ) ; // Structured spatial effects (CAR)
	PARAMETER(           log_tau ) ; // Log-precision for random effects. Passed as log(tau) so that tau is enforced positive.
	PARAMETER(       logodds_phi ) ; // Log-odds of mixture proportion of spatial vs independent components.
	PARAMETER(            log_nu ) ; // nu parameter of generalised logit, encoded as log to keep positive

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
	Type nu = exp( log_nu ) ;

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
	// or a an exponential prior on the sd.
	// We found dexp() on the sd led to numerical errors so this is a
	// direct implementation of the Type 2 Gumbel distribution:
	//nll -= log_dN( log_tau, 0.0, prior_sd_rate ) ;
	nll -= (
		log(prior_sd_rate) - log(2.0) - (3.0/2.0)*log_tau - (prior_sd_rate / sqrt(tau))
	) ;
	//nll -= dexp( sd_of_random_effects, prior_sd_rate, /* give_log */ 1 ) ;

	// Normal prior on logodds phi
	nll -= log_dN( logodds_phi, prior_logodds_phi_mean, prior_logodds_phi_sd ) ;

	//dexp( sd_of_random_effects, prior_sd_rate, /* give_log */ 1 ) ;
	nll -= ( log_dN( u, 0.0, 1.0 ).array() * N.array() ).sum() ; // independent random effects are standard gaussian

	// Prior on nu.  nu = 1 is plain logistic so put prior centred on this.
	nll -= log_dN( log_nu, 0.0, prior_log_nu_sd ) ;

	// ICAR (spatial) penalty term
	vector<Type> Qv = Q * v;
	Type penalty = -Type(0.5) * (v * Qv).sum() ;
	nll -= penalty ;

	// Linear predictor and binomial likelihood
	vector<Type> predictor = intercept + (beta * x) ;
	if( model_choice == "bym2" ) {
		predictor += (
			sd_of_random_effects * sqrt(1-phi)*u
			+ sd_of_random_effects * sqrt(phi)*v
		) ;
	}
	vector<Type> p ;
	if( link_choice == "logit" ) {
		p = invlogit( predictor ) ;
	} else if( link_choice == "generalised-logit" ) {
		// logit with additional nu parameter
		// https://en.wikipedia.org/wiki/Generalised_logistic_function
		p = (
			Type(1.0) / pow( (Type(1.0) + exp( -predictor )), (Type(1.0)/nu) )
		) ;
	} else if( link_choice == "linear" ) {
		p = clamp( predictor, 0.001, 0.999 ) ;
	}
	nll -= sum( dbinom(y, N, p, true) ); 

	return nll;
}
