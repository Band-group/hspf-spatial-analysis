suppressMessages( library( dplyr ))
suppressMessages( library( argparse ))
suppressMessages( library( stringr ))
suppressMessages( library( readr ))
suppressMessages( library( tibble ))
suppressMessages( library( tools))

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
		help = "grid cell size. FIXME: would be better to include in the fit object!",
		required = TRUE
	)
	parser$add_argument(
		"--output",
		type = "character",
		help = "Filename of .tsv file to write"
	)
	parser$add_argument(
		"--hspf_covariates",
		type = "character",
		help = "covariates used in hspf model",
		required = TRUE,
		default = "pfpr2000"
	)
	return( parser$parse_args() )
}

options( width = 300 )
args = parse_arguments()
#testing
# args <- list()
# args$area = "DRC"
# args$fit = c("output/pf=pf8-version/hspf/fixed-r0=25.0-sigma0=0.6-fc=none/grid-type=hexagon-size=1/Pfsa3/Pfsa3-model=bym2+fc=none-200km-area=DRC-min_N=0.rds")
# args$output = "output/pf=pf8-version/hspf/fixed-r0=25.0-sigma0=0.6-fc=none/grid-type=hexagon-size=1/Pfsa3/Pfsa3-model=bym2+fc=none-200km-area=DRC-min_N=0-summary.tsv"
# args$min_N = 0
# args$cellsize = 1
# args$hspf_covariates = "none"


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
		#cellsize = fit$cellsize, #FIXME
		cellsize = args$cellsize,
		HbSr0 = fit$r0,
		HbSsigma0 = fit$sigma0,
		allele = fit$allele,
		area = args$area,
		countries = paste( fit$areas, collapse = "," ),
		min_km_to_survey_pt = fit$min_km_to_survey_pt,
		min_N = args$min_N,
		covariate = args$hspf_covariates,
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

# create/overwrite
result$Reported <- NA_character_

# move to first column
result <- result[, c("Reported", setdiff(names(result), "Reported")), drop = FALSE]

result$Reported <- ifelse(
  result$celltype == "hexagon" &
    result$cellsize == 1 &
    result$HbSr0 == 25 &
    result$HbSsigma0 == 0.6 &
    result$area %in% c("africa") &
    result$allele %in% c("Pfsa1","Pfsa3") &
    result$covariate == "pfpr2000",
  "Figure 1, Figure 2",
  ifelse(
    result$celltype == "hexagon" &
      result$cellsize == 1 &
      result$HbSr0 == 25 &
      result$HbSsigma0 == 0.6 &
      result$area %in% c("mauritania","senegal+gambia","ghana","nigeria","drc","DRC","uganda","tanzania","mozambique") &
      result$allele %in% c("Pfsa1") &
      result$covariate == "pfpr2000",
    "Figure S3",
    ifelse(
      result$celltype == "hexagon" &
        result$cellsize == 1 &
        result$HbSr0 == 25 &
        result$HbSsigma0 == 0.6 &
        result$area %in% c("global", "DRC+east","africa","drc+east") &
        result$allele %in% c("Pfsa1","Pfsa2","Pfsa13","Pfsa4") &
        result$covariate == "pfpr2000",
      "Figure 2",
      "Table S3 only"
    )
  )
)
result <- result %>%
  mutate(
    area = recode(
      area,
      "waf" = "West Africa",
      "drc+east" = "Central and Eastern Africa",
      "eaf" = "East Africa",
	  "gambia+senegal" = "Gambia and Senegal",
	  "global" = "Global"
    ),
    area = tools::toTitleCase(area)  # Capitalizes first letter of each word
  )
echo( "++ Ok, saving results as a .tsv table to %s...\n", args$output )
readr::write_tsv( result, args$output, append = file.exists( args$output ))

echo( "++ Success.\n" )
echo( "++ Thank you for using summarise_hspf_fits.R!\n" )
