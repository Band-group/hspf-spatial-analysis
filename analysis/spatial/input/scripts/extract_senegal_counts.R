library( dplyr )

library( argparse )

echo = function( message, ... ) {
	cat( sprintf( message, ... ))
}

parse_arguments <- function() {
	parser = ArgumentParser(
		description = 'Extract Pfsa counts'
	)
	parser$add_argument(
		"--indir",
		type = "character",
		help = "path to folder containing Pf7 data",
		default = "../../../data/senegal"
	)
	parser$add_argument(
		"--output",
		type = "character",
		help = "path to output directory",
		default = "input/hbs-pf-v2.sqlite",
		required = TRUE
	)
	
	return( parser$parse_args() )
}

args = parse_arguments()

paths = list(
	dosage = sprintf( "%s/calls.dosage.gz", args$indir ),
	sampmap = sprintf( "%s/all_sampmap.txt", args$indir ),
	sites = sprintf( "%s/sites.tsv", args$indir )
)

#functions
load.entry.from.Rdata <- function( filename, what ) {
  env = new.env()
  load( file = filename, envir = env )
  # Sanity check - we need these:
  stopifnot( what %in% names(env))
  result = env[[what]]
  rm(env)
  return( result )
}

dosage = readr::read_table( paths$dosage )
sample_map = readr::read_tsv( paths$sampmap )
#samples = readr::read_csv( paths$samples )
#annotated = readr::read_tsv( paths$annotated )
#dupes = readr::read_table( paths$dupes, col_names = "sample" )
sites = readr::read_tsv( paths$sites, comment = '#' )

variants = dosage[,1:6]
genotypes = as.matrix(dosage[,7:ncol(dosage)])
# Remove any mixed genotypes
genotypes[ genotypes == 1 ] = NA

# Fix missing allele calls
variants$alleleB[ variants$position == 814288 ] = "T"
variants$alleleB[ variants$position == 814329 ] = "G"

rownames( genotypes ) = sprintf( "%s:%d:%s>%s", variants$chromosome, variants$position, variants$alleleA, variants$alleleB )

result = tibble::tibble(
	old_sample_name = colnames( genotypes ),
	new_sample_name = sample_map$new_samp_name[ match( colnames(genotypes), sample_map$old_samp_name )]
)
result$country = stringr::str_sub( result$new_sample_name, 1, 3 )
result$site = stringr::str_sub( result$new_sample_name, 5, 7 )
result$year = stringr::str_sub( result$new_sample_name, 9, 12 )
result = (
	result
	%>% inner_join( sites %>% select( site = Site, longitude, latitude ), by = "site" )
)
int = as.integer
echo( "++ Computing by_sample..." )
by_sample = dplyr::bind_cols(
	result,
	t(genotypes)
) %>% mutate(
	source = "Schaffner et al Senegal 2023",
	study = "https://doi.org/10.1038/s41467-023-43087-4",
	datatype = 'WGS',
	country = "Senegal",
	year,
	N = 1,
	`Pfsa1:ref` = int(`Pf3D7_02_v3:631190:T>A` == 0),
	`Pfsa1:nonref` = int(`Pf3D7_02_v3:631190:T>A` == 2),
	`Pfsa2:ref` = int(`Pf3D7_02_v3:814288:C>T` == 0),
	`Pfsa2:nonref` = int(`Pf3D7_02_v3:814288:C>T` == 2),
	`Pfsa3:ref` = int(`Pf3D7_11_v3:1058035:T>A` == 0 ),
	`Pfsa3:nonref` = int(`Pf3D7_11_v3:1058035:T>A` == 2 ),
	`Pfsa4:ref` = int(`Pf3D7_04_v3:1121472:T>A` == 0),
	`Pfsa4:nonref` = int(`Pf3D7_04_v3:1121472:T>A` == 2 ),
	exclude = 'no'
) %>% select(
	source,
	study,
	datatype,
	country,
	year,
	site,
	latitude,
	longitude,
	ID = new_sample_name,
	N,
	`Pfsa1:ref`,
	`Pfsa1:nonref`,
	`Pfsa2:ref`,
	`Pfsa2:nonref`,
	`Pfsa3:ref`,
	`Pfsa3:nonref`,
	`Pfsa4:ref`,
	`Pfsa4:nonref`,
	exclude
)

echo( "++ Computing by_site..." )
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

echo( "++ Outputting...\n" )
db = DBI::dbConnect( RSQLite::SQLite(), args$output )
tables = DBI::dbGetQuery( db, "SELECT * FROM sqlite_master WHERE type == 'table'" )
print( tables )
if( "by_site" %in% tables$name ) {
	DBI::dbExecute( db, "DELETE FROM by_site WHERE source == 'schaffner_et_al_2023' ")
}
DBI::dbWriteTable( db, "by_site", by_site, overwrite = FALSE, append = TRUE )
if( "by_sample" %in% tables$name ) {
	DBI::dbExecute( db, "DELETE FROM by_sample WHERE source == 'schaffner_et_al_2023' ")
}
DBI::dbWriteTable( db, "by_sample", by_sample, overwrite = FALSE, append = TRUE )
DBI::dbDisconnect( db )
