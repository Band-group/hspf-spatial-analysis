library( tidyverse )
library( dplyr )
library( dbplyr )
library( RSQLite )
library( argparse )

parse_arguments <- function() {
	parser = ArgumentParser(
		description = 'Extract Pfsa counts'
	)
	parser$add_argument(
		"--indir",
		type = "character",
		help = "path to folder containing DRC input data",
		default = "input/dr_congo"
	)
	parser$add_argument(
		"--output",
		type = "character",
		help = "path to output directory",
		default = "input/hbs-pf-v4.sqlite",
		required = TRUE
	)
	
	return( parser$parse_args() )
}

args = parse_arguments()

paths = list(
	data = sprintf( "%s/biallelic_processed0.rds", args$indir )
)

data = readRDS( paths$data )
stopifnot( length( which( rownames( data$samples ) != rownames( data$counts ))) == 0 )
stopifnot( length( which( rownames( data$samples ) != rownames( data$coverage ))) == 0 )

SNPs = tibble(
	chromosome = sprintf( "chr%d", c( 2, 2, 2, 4, 4, 11, 11 )),
	position = c(
		629996, 631190,
		814329,
		1121472, 1122147,
		1058035, 1057437
	),
	locus = c(
		rep( "Pfsa1", 2 ),
		rep( "Pfsa2", 1 ),
		rep( "Pfsa4", 2 ),
		rep( "Pfsa3", 2 )
	),
	type = c(
		"secondary", "lead",
		"lead",
		"lead", "secondary",
		"lead", "lead"
	)
)

samples = data$samples
coverage = data$coverage[, sprintf( "%s_%d", SNPs$chromosome, SNPs$position )]
counts = data$counts[, sprintf( "%s_%d", SNPs$chromosome, SNPs$position )]

calls = matrix(
	nrow = nrow( samples ),
	ncol = ncol( counts ),
	dimnames = list(
		rownames( counts ),
		colnames( counts )
	)
)

threshold = 0.9
calls[,] = NA
calls[ counts / coverage >= threshold ] = 1
calls[ counts / coverage <= (1-threshold) ] = 0
calls[ coverage < 5 ] = NA
table( calls[,2], is.na( coverage[,2] ), useNA="always" )

ratios = matrix(
	nrow = nrow( samples ),
	ncol = ncol( counts ),
	dimnames = list(
		rownames( counts ),
		sprintf( "%s_ratio", colnames( counts ))
	)
)
ratios[,] = counts/coverage

samples = as_tibble( bind_cols( samples, calls, ratios ))
samples = (
	samples
	%>% mutate(
		latitude = lat,
		longitude = long,
		source = "Verity et al 2021",
		study = "Verity et al 2021",
		datatype = "MIP",
		country = c(
			'DRC' = 'Democratic_Republic_of_the_Congo',
			'Ghana' = 'Ghana',
			'Tanzania' = 'Tanzania',
			'Uganda' = 'Uganda',
			'Zambia' = 'Zambia'
		)[Country],
		year = Year,
		site = NA,
		N = 1,
		# We have figured out the alleles are encoded the other way round.
		# c.f. email trail with Bob Verity 12th Feb 2025
		# Fixing this here.
		`Pfsa1:ref` = `chr2_631190`,
		`Pfsa1:nonref` = 1 - `chr2_631190`,
		`Pfsa2:ref` = `chr2_814329`,
		`Pfsa2:nonref` = 1 - `chr2_814329`,
		`Pfsa3:ref` = `chr11_1058035`,
		`Pfsa3:nonref` = 1 - `chr11_1058035`,
		`Pfsa4:ref` = `chr4_1121472`,
		`Pfsa4:nonref` = 1 - `chr4_1121472`,
		exclude = "no"
	)
)
# If fklipped alleles fixes this, don't need to exclude
#samples$exclude[ samples$country != 'Democratic_Republic_of_the_Congo' ] = 'yes'

by_sample = (
	samples
	%>% select(
		source, study, datatype, country, year, site, latitude, longitude,
		ID, N,
		`Pfsa1:ref`, `Pfsa1:nonref`,
		`Pfsa2:ref`, `Pfsa2:nonref`,
		`Pfsa3:ref`, `Pfsa3:nonref`,
		`Pfsa4:ref`, `Pfsa4:nonref`,
		`exclude`
	)
)

by_site = (
	by_sample
	%>% filter( !is.na( latitude ))
	%>% filter( exclude == 'no' )
	%>% group_by(
		source, study, datatype, country, year, site, latitude, longitude
	) %>% summarise(
		N = n(),
		'Pfsa1:ref' = sum(`Pfsa1:ref`, na.rm = T ), `Pfsa1:nonref` = sum( `Pfsa1:nonref`, na.rm = T ),
		'Pfsa2:ref' = sum(`Pfsa2:ref`, na.rm = T ), `Pfsa2:nonref` = sum( `Pfsa2:nonref`, na.rm = T ),
		'Pfsa3:ref' = sum(`Pfsa3:ref`, na.rm = T ), `Pfsa3:nonref` = sum( `Pfsa3:nonref`, na.rm = T ),
		'Pfsa4:ref' = sum(`Pfsa4:ref`, na.rm = T ), `Pfsa4:nonref` = sum( `Pfsa4:nonref`, na.rm = T ),
		exclude = max( exclude )
	)
)

options(width=200)
print( by_site )

db = DBI::dbConnect( RSQLite::SQLite(), args$output )
DBI::dbExecute( db, "DELETE FROM by_site WHERE source == 'Verity et al 2021' ")
DBI::dbWriteTable( db, "by_site", by_site, overwrite = FALSE, append = TRUE )
DBI::dbExecute( db, "DELETE FROM by_sample WHERE source == 'Verity et al 2021' ")
DBI::dbWriteTable( db, "by_sample", by_sample, overwrite = FALSE, append = TRUE )
DBI::dbDisconnect( db )
