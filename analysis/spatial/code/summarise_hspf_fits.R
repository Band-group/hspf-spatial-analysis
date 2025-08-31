suppressMessages( library( dplyr ))
suppressMessages( library( argparse ))
suppressMessages( library( stringr ))

source( "code/figures/fig1_impl.R" )

echo <- function( message, ... ) {
	cat( sprintf( message, ... ))
}

parse_arguments <- function() {
	parser = ArgumentParser(
		description = 'Collate hspf fit output into one file'
	)
	parser$add_argument(
		"--fit",
		type = "character",
		nargs = "+",
		help = "Fit object, from BYM.R, to summarise",
		required = TRUE
	)
	parser$add_argument(
		"--area",
		type = "character",
		help = "Area specification",
		required = TRUE
	)
	parser$add_argument(
		"--min_N",
		type = "character",
		help = "min N. FIXME: would be better to include in the fit object!",
		required = TRUE
	)
	parser$add_argument(
		"--cellsize",
		type = "character",
		help = "grid cell size.  FIXME: would be better to include in the fit object!",
		required = TRUE
	)
	parser$add_argument(
		"--output",
		type = "character",
		help = "Filename of .tsv file to write"
	)
	return( parser$parse_args() )
}

options( width = 300 )
args = parse_arguments()

for( filename in args$fit ) {
	if( !file.exists( filename )) {
		stop( "!! File %s not found, quitting.\n", filename )
	}
}

result = tibble()
echo( "  ... processing %s...\n", args$fit )
fit = readRDS( args$fit )
fit$sampled_parameters$posterior.sample = 1:nrow( fit$sampled.parameters )

echo( "++ Predicting...\n" )
predictions = make_hspf_curves(
	fit$sampled.parameters,
	at = c( 0.1, 0.2 ),
	link_fn = list(
		logit = function( v, parameters ) {
			x = parameters[['intercept']] + parameters[['beta']]*v
			return( exp(x)/(1+exp(x)) )
		},
		`generalised-logit` = function( v, parameters ) {
			x = parameters[['intercept']] + parameters[['beta']]*v
			nu = exp( parameters[['log_nu']] )
			return( 1/(1 + exp(-x))^(1/nu))
		},
		linear = function( v, parameters ) {
			x = parameters[['intercept']] + parameters[['beta']]*v
			return( pmax( pmin( x, 0.999 ), 0.001 ))
		}
	)[[fit$link]]
)
compute.delta = function( x, y ) {
	return( y[x == 0.2] - y[x == 0.1 ])
}
delta_summary = (
	predictions
	%>% group_by(
		posterior.sample
	)
	%>% summarise(
		delta = compute.delta( x, y )
	)
	%>% ungroup()
	%>% summarise(
		delta_mean = mean(delta),
		delta_median = median(delta),
		delta_q2.5 = quantile( delta, p = 0.025 ),
		delta_q97.5 = quantile( delta, p = 0.975 )
	)
)

echo( "++ Summarising...\n" )
summary = bind_cols(
	delta_summary,
	fit$sampled.parameters
	%>% summarise(
		pf_at_0.05 = mean( gl( 0.05, pick( intercept, beta, log_nu)), na.rm = T ),
		pf_at_0.1 = mean( gl( 0.1, pick( intercept, beta, log_nu)), na.rm = T ),
		pf_at_0.1.q2.5 = quantile( gl( 0.1, pick( intercept, beta, log_nu)), 0.025 ),
		pf_at_0.1.q97.5 = quantile( gl( 0.1, pick( intercept, beta, log_nu)), 0.975 ),
		pf_at_0.15 = mean( gl( 0.15, pick( intercept, beta, log_nu)), na.rm = T ),
		pf_at_0.2 = mean( gl( 0.2, pick( intercept, beta, log_nu)), na.rm = T ),
		pf_at_0.2.q2.5 = quantile( gl( 0.2, pick( intercept, beta, log_nu)), 0.025 ),
		pf_at_0.2.q97.5 = quantile( gl( 0.2, pick( intercept, beta, log_nu)), 0.975 ),
		pf_at_0.25 = mean( gl( 0.25, pick( intercept, beta, log_nu)), na.rm = T ),
		pf_at_0.3 = mean( gl( 0.3, pick( intercept, beta, log_nu)), na.rm = T ),
		beta.mean = mean( beta, na.rm = TRUE ),
		beta.q2.5 = quantile( beta, 0.025 ),
		beta.q25 = quantile( beta, 0.25 ),
		beta.q50 = quantile( beta, 0.5 ),
		beta.q75 = quantile( beta, 0.75 ),
		beta.q97.5 = quantile( beta, 0.975 )
	)
)

echo( "++ Forming result...\n" )
print( paste( fit$areas, collapse = "," ) )
result = bind_cols(
	tibble(
		celltype = fit$celltype,
#		cellsize = fit$cellsize, # FIXME
		cellsize = args$cellsize,
		HbSr0 = fit$r0,
		HbSsigma0 = fit$sigma0,
		allele = fit$allele,
		area = args$area,
		countries = paste( fit$areas, collapse = "," ),
		min_km_to_survey_pt = fit$min_km_to_survey_pt,
		min_N = args$min_N,
		model = fit$model,
		transform = fit$transform,
		n_data_points = nrow( fit$data ),
		mean_cpo = mean( fit$summary$cpo,na.rm=TRUE ),
		mean_waic = mean( fit$summary$waic,na.rm=TRUE ),
		mean_ll_integrated = mean( fit$summary$marginal_ll_integration,na.rm=TRUE ),
		mean_ll_gaussian = mean( fit$summary$marginal_ll_gaussian,na.rm=TRUE )
	),
	summary
)

echo( "++ Ok, saving results as a .tsv table to %s...\n", args$output )
readr::write_tsv( result, args$output, append = file.exists( args$output ))

#echo( "++ Re-reading results from %s...\n", args$output )

#updatedresult <- readr::read_tsv(args$output)

use_tex = FALSE
if( use_tex ) {
	suppressMessages( library( knitr ))
	suppressMessages( library( kableExtra ))

	echo( "++ Ok, %d rows read.\n", nrow( updatedresult ))
	# Save table as latex
	resulttex <- updatedresult %>%
	select(-countries) %>%
	mutate(across(where(is.numeric), ~ round(., 2)))
	# Remove underscore in pf allele (generate issues later)
	resulttex$allele <- gsub("_", "", resulttex$allele)
	resulttex$allele <- gsub("\\+", "", resulttex$allele)

	#get unique range and sigma for file name 
	addsuffix <- paste0("r0=",unique(fit$r0),"-sigma0=", unique(fit$sigma0))

	# Rename columns for better readability
	colnames(resulttex) <- c("Type", "Size", "$\\rho_0$","$\\sigma_0$","Pf allele", 
						"Domain", "km" ,"Model", "Link", "N", "CPO", "WAIC", "LL(Int.)", 
						"LL(Gauss.)", "Mean", 
						"Q2.5", "Q25", "Median", "Q75", "Q97.5")
	#remove some variables
	resulttex <- resulttex %>%
	select(-c(Q25,Q75))
	#rename some entries
	resulttex <- resulttex %>% 
	mutate(across('Model', stringr::str_replace, 'norandom', 'Besag')) %>% 
	mutate(across('Model', stringr::str_replace, 'bym2', 'BYM')) %>% 
	mutate(across('Domain', stringr::str_replace, 'waf', 'West Af.')) %>% 
	mutate(across('Domain', stringr::str_replace, 'eaf', 'East Af.'))

	resulttex$Domain <- resulttex$Domain %>% stringr::str_to_title()
	resulttex$Link <- resulttex$Link %>% stringr::str_to_title()
	resulttex$Type <- resulttex$Type %>% stringr::str_to_title()

	# Create a formatted table with multicolumns
	captiontext <- "Results summary. Assessing the effects of HbS allele frequency on P.falciparum (Pf) allele frequency in Africa
	(West Africa: Gambia, Senegal, Mali, Benin, Burkina Faso, Ivory Coast, Ghana, Guinea, Mauritania, Nigeria, Senegal, Togo, 
	Angola, Cameroon; East Africa: Ethiopia, Kenya, Madagascar, Malawi, Mozambique, Rwanda, Uganda, and United Republic of Tanzania; 
	Gabon; Other African countries: Central African Republic, Republic of the Congo, Democratic Republic of the Congo). 
	Models are specified as follows: Spatial unit (cell) type and size (hexagon: distance between oppositve edges; square: edge length. Unit: degree) of the discretised spatial domain;
	Penalised complexity prior thresholds on the spatial random field (range and sigma) for the HbS models; 
	Pf allele locus; Study domain; Radius distance used to subset spatial units; Type of spatial Pf model (Besag includes only structured effects; BYM includes both structured and unstructured effects); Pf model link function; 
	Number of observations in Pf model; Predictive performance metrics including Conditional predictive ordinates (CPO), Watanabe-Akaike information criterion (WAIC), log likelihood 
	estimated with integrative or Gaussian approximations; Estimated values (mean, median, quantiles 0.025, 0.5 (median), 0.975) of the HbS effects on Pf allele frequency."

	echo( "++ Ok, saving results as a .tex table to %s...\n", args$output )
	resulttex <- resulttex %>%
	kable("latex", booktabs = TRUE, caption = captiontext, longtable = TRUE,
			escape=FALSE,linesep=c("", "", "","", "", "","", "\\addlinespace")) %>%
	kable_classic(full_width = F) %>%
	add_header_above(c("Spatial unit" = 2, "HbS priors" = 2, rep("",6), "Predictive Performance" = 4, "Estimated HbS effect" = 4)) %>%
	kable_styling(latex_options = c("repeat_header"), font_size = 7,
					repeat_header_continued = "\\textit{(Continued on next page...)}") 
	# Save the table to a .tex file
	resulttex_name = gsub( "[.]tsv$", "", args$output)
	save_kable(resulttex, file = paste0(resulttex_name,"_",addsuffix,".tex"))
}

echo( "++ Success.\n" )
echo( "++ Thank you for using summarise_hspf_fits.R!\n" )
