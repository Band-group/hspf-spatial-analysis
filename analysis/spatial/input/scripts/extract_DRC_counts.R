library( tidyverse )
library( dplyr )
library( dbplyr )
library( rbgen )

data = readRDS( "data/dr_congo/biallelic_processed0.rds" )
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

calls[ counts / coverage >= 0.9 ] = 1
calls[ counts / coverage <= 0.1 ] = 0
calls[ coverage < 5 ] = NA
table( calls[,2] )

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
samples$latitude = samples$lat
samples$longitude = samples$long

by_site = (
	samples
	%>% group_by(
		STUDY_CODE, Country, latitude, longitude
	) %>% summarise(
		N = n(),
		`Pfsa1:ref` = sum( !is.na( `chr2_631190` )) - sum( `chr2_631190`, na.rm = T ),
		`Pfsa1:nonref` = sum( `chr2_631190`, na.rm = T ),
		`Pfsa1:ref2` = sum( !is.na( `chr2_631190_ratio` )) - sum( `chr2_631190_ratio`, na.rm = T ),
		`Pfsa1:nonref2` = sum( `chr2_631190_ratio`, na.rm = T ),
		`Pfsa2:ref` = sum( !is.na( `chr2_814329` )) - sum( `chr2_814329`, na.rm = T ),
		`Pfsa2:nonref` = sum( `chr2_814329`, na.rm = T ),
		`Pfsa2:ref2` = sum( !is.na( `chr2_814329_ratio` )) - sum( `chr2_814329_ratio`, na.rm = T ),
		`Pfsa2:nonref2` = sum( `chr2_814329_ratio`, na.rm = T ),
		`Pfsa3:ref` = sum( !is.na( `chr11_1058035` )) - sum( `chr11_1058035`, na.rm = T ),
		`Pfsa3:nonref` = sum( `chr11_1058035`, na.rm = T ),
		`Pfsa3:ref2` = sum( !is.na( `chr11_1058035_ratio` )) - sum( `chr11_1058035_ratio`, na.rm = T ),
		`Pfsa3:nonref2` = sum( `chr11_1058035_ratio`, na.rm = T ),
		`Pfsa4:ref` = sum( !is.na( `chr4_1121472` )) - sum( `chr4_1121472`, na.rm = T ),
		`Pfsa4:nonref` = sum( `chr4_1121472`, na.rm = T ),
		`Pfsa4:ref2` = sum( !is.na( `chr4_1121472_ratio` )) - sum( `chr4_1121472_ratio`, na.rm = T ),
		`Pfsa4:nonref2` = sum( `chr4_1121472_ratio`, na.rm = T )
	)
)
by_site$source = "Verity_et_al_2021"
by_site$study = "Verity_et_al_2021"
by_site$country = c(
	'DRC' = 'Democratic_Republic_of_the_Congo',
	'Ghana' = 'Ghana',
	'Tanzania' = 'Tanzania',
	'Uganda' = 'Uganda',
	'Zambia' = 'Zambia'
)[by_site$Country]
by_site$site = NA
by_site = by_site[,
	c(
		"source", "study", "country", "site", "latitude", "longitude",
		"N",
		"Pfsa1:ref", "Pfsa1:nonref",
		"Pfsa1:ref2", "Pfsa1:nonref2",
		"Pfsa2:ref", "Pfsa2:nonref",
		"Pfsa2:ref2", "Pfsa2:nonref2",
		"Pfsa3:ref", "Pfsa3:nonref",
		"Pfsa3:ref2", "Pfsa3:nonref2",
		"Pfsa4:ref", "Pfsa4:nonref",
		"Pfsa4:ref2", "Pfsa4:nonref2"
	)
]

write_tsv( samples, "results/genotypes/Verity_et_al_2021_by_sample.tsv" )
write_tsv( by_site, "results/genotypes/Verity_et_al_2021_by_site.tsv" )

db = DBI::dbConnect( RSQLite::SQLite(), "results/genotypes/hbs-pf.sqlite" )
DBI::dbExecute( db, "DELETE FROM by_site WHERE source == 'Verity_et_al_2021' ")
DBI::dbWriteTable( db, "by_site", by_site, overwrite = FALSE, append = TRUE )
DBI::dbDisconnect( db )
