library( rnaturalearth )
library( rnaturalearthhires )
library( dplyr )

paths = list(
	dosage = "data/senegal/reads/fetched/calls.variants.dosage.gz",
	sampmap = "data/senegal/all_sampmap.txt",
	dupes = "data/senegal/all_sampdupes.txt",
	annotated = "data/senegal/PRJNA972644_annotated.tsv",
	samples = "data/senegal/PRJNA972644.csv",
	sites = "data/senegal/sites.tsv",
	senegal = "data/senegal/shapes/sen_admbnda_adm2_anat_20240520.shp"
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
samples = readr::read_csv( paths$samples )
annotated = readr::read_tsv( paths$annotated )
dupes = readr::read_table( paths$dupes, col_names = "sample" )
senegal = sf::read_sf( paths$senegal )
sites = readr::read_tsv( paths$sites, comment = '#' )

variants = dosage[,1:6]
genotypes = dosage[,7:ncol(dosage)]
rownames( genotypes ) = sprintf( "%s:%d:%s>%s", variants$chromosome, variants$position, variants$alleleA, variants$alleleB )

result = tibble::tibble(
	old_sample_name = colnames( genotypes ),
	new_sample_name = sample_map$new_samp_name[ match( colnames(genotypes), sample_map$old_samp_name )]
)
result$country = stringr::str_sub( result$new_sample_name, 1, 3 )
result$site = stringr::str_sub( result$new_sample_name, 5, 7 )
result$year = stringr::str_sub( result$new_sample_name, 9, 12 )
sites = readr::read_tsv( paths$sites, comment = '#' )
result = (
	result
	%>% inner_join( sites %>% select( site = Site, longitude, latitude ), by = "site" )
)
int = as.integer
by_sample = dplyr::bind_cols(
	result,
	t(genotypes)
) %>% mutate(
	source = "schaffner_et_al_2023",
	study = "schaffner_et_al_2023",
	country = "Senegal",
	N = 1,
	`Pfsa1:ref` = int(`Pf3D7_02_v3:631190:T>A` == 0),
	`Pfsa1:nonref` = int(`Pf3D7_02_v3:631190:T>A` == 2),
	`Pfsa2:ref` = NA,
	`Pfsa2:nonref` = NA,
	`Pfsa3:ref` = int(`Pf3D7_11_v3:1058035:T>A` == 0 ),
	`Pfsa3:nonref` = int(`Pf3D7_11_v3:1058035:T>A` == 2 ),
	`Pfsa4:ref` = int(`Pf3D7_04_v3:1121472:T>A` == 0),
	`Pfsa4:nonref` = int(`Pf3D7_04_v3:1121472:T>A` == 2 ),
	exclude = 'no'
) %>% select(
	source,
	study,
	country,
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

by_site = (
	by_sample
	%>% group_by(
		source, study, country, site, latitude, longitude
	) %>% summarise(
		N = n(),
		'Pfsa1:ref' = sum(`Pfsa1:ref`, na.rm = T ), `Pfsa1:nonref` = sum( `Pfsa1:nonref`, na.rm = T ),
		'Pfsa2:ref' = sum(`Pfsa2:ref`, na.rm = T ), `Pfsa2:nonref` = sum( `Pfsa2:nonref`, na.rm = T ),
		'Pfsa3:ref' = sum(`Pfsa3:ref`, na.rm = T ), `Pfsa3:nonref` = sum( `Pfsa3:nonref`, na.rm = T ),
		'Pfsa4:ref' = sum(`Pfsa4:ref`, na.rm = T ), `Pfsa4:nonref` = sum( `Pfsa4:nonref`, na.rm = T ),
		exclude = max( exclude )
	)
)

db = DBI::dbConnect( RSQLite::SQLite(), "results/genotypes/hbs-pf-v2.sqlite" )
DBI::dbExecute( db, "DELETE FROM by_site WHERE source == 'schaffner_et_al_2023' ")
DBI::dbWriteTable( db, "by_site", by_site, overwrite = FALSE, append = TRUE )
DBI::dbExecute( db, "DELETE FROM by_sample WHERE source == 'schaffner_et_al_2023' ")
DBI::dbWriteTable( db, "by_sample", by_sample, overwrite = FALSE, append = TRUE )
DBI::dbDisconnect( db )
