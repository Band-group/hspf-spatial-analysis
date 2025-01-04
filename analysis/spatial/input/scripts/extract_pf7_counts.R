library( tidyverse )
library( dplyr )
library( dbplyr )
library( rbgen )

library( argparse )

options(width=200)
echo <- function( message, ... ) {
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
		default = "/well/band/projects/pf7/"
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
	samples = sprintf( "%s/data/samples/Pf7_samples.txt", args$indir ),
	genotypes = sprintf( "%s/results/bgen/pf7.filtered.bgen", args$indir )
)

load.genotypes <- function( filename, SNPs ) {
	G = bgen.load(
		filename,
		ranges = tibble(
			chromosome = SNPs$chromosome,
			start = SNPs$position,
			end = SNPs$position
		),
		max_entries = 10
	)
	# for the sake of this analysis, we ignore mixed calls
	# or calls of other allels
	G$data = G$data[,,c(1,3)]
	G$dosage = G$data[,,2]
	G$dosage[ rowSums(G$data,dims=2) == 0 ] = NA

	G$variants$allele0[1] = "C"
	G$variants$allele1[1] = "T"
	G$variants$allele0[9] = "T"
	G$variants$allele1[9] = "A"
	G$variants$ID = sprintf(
		"chr%d:%d:%s>%s",
		as.integer( gsub( "_v3", "", gsub( "Pf3D7_", "", G$variants$chromosome ))),
		G$variants$position,
		G$variants$allele0,
		G$variants$allele1
	)
	return(G)
}

samples = read_tsv( paths$samples )
samples$Country = gsub( " ", "_", samples$Country )
samples$Country[ grep( "Ivoire", samples$Country)] = "Cote_dIvoire" 

SNPs = tibble(
	chromosome = sprintf( "Pf3D7_%02d_v3", c( 2, 2, 2, 2, 2, 4, 4, 11, 11 )),
	position = c(
		629996, 631190, 630290,
		814288, 814329,
		1121472, 1122147,
		1058035, 1057437
	),
	locus = c(
		rep( "Pfsa1", 3 ),
		rep( "Pfsa2", 2 ),
		rep( "Pfsa4", 2 ),
		rep( "Pfsa3", 2 )
	),
	type = c(
		"secondary", "lead", "secondary",
		"lead", "secondary",
		"lead", "secondary",
		"lead", "lead"
	)
)

G = load.genotypes( paths$genotypes, SNPs )

for( i in 1:nrow( G$variants )) {
	samples[,G$variants$ID[i]] = G$dosage[i,]
}

wIn = which( samples$`Exclusion reason` == 'Analysis_set' )
samples$exclude = 'yes'
samples$exclude[wIn] = 'no'
#write_tsv( samples, "results/pfsa_genotypes/pf7.tsv" )

samples$source = "MalariaGEN Pf7"

samples = (
	samples
	%>% mutate(
		ID = Sample,
		latitude =  `Admin level 1 latitude`,
		longitude = `Admin level 1 longitude`,
		source = "MalariaGEN Pf7",
		study = Study,
		datatype = "WGS",
		country = Country,
		site = `Admin level 1`,
		`Pfsa1:ref` = 1 - `chr2:631190:T>A`,
		`Pfsa1:nonref` = `chr2:631190:T>A`,
		`Pfsa2:ref` = 1 - `chr2:814288:C>T`,
		`Pfsa2:nonref` = `chr2:814288:C>T`,
		`Pfsa3:ref` = 1 - `chr11:1058035:T>A`,
		`Pfsa3:nonref` = `chr11:1058035:T>A`,
		`Pfsa4:ref` = 1 - `chr4:1121472:T>A`,
		`Pfsa4:nonref` = `chr4:1121472:T>A`,
		exclude = 'no'
	)
)

by_sample = (
	samples
	%>% mutate(
		ID = Sample,
		N = 1,
	)
	%>% select(
		source, study, datatype, country, site, latitude, longitude,
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
		source, study, datatype, country, site, latitude, longitude
	)
	%>% summarise(
		N = sum(N),
		`Pfsa1:ref` = sum( `Pfsa1:ref`, na.rm = T ), `Pfsa1:nonref` = sum( `Pfsa1:nonref`, na.rm = T ),
		`Pfsa2:ref` = sum( `Pfsa2:ref`, na.rm = T ), `Pfsa2:nonref` = sum( `Pfsa2:nonref`, na.rm = T ),
		`Pfsa3:ref` = sum( `Pfsa3:ref`, na.rm = T ), `Pfsa3:nonref` = sum( `Pfsa3:nonref`, na.rm = T ),
		`Pfsa4:ref` = sum( `Pfsa4:ref`, na.rm = T ), `Pfsa4:nonref` = sum( `Pfsa4:nonref`, na.rm = T ),
		`exclude` = 'no'
	)
)

print( by_site )

echo( "++Outputting to %s...\n", args$output )
	db = DBI::dbConnect( RSQLite::SQLite(), args$output )
	DBI::dbExecute( db, "DELETE FROM by_site WHERE source == 'MalariaGEN Pf7' ")
	DBI::dbWriteTable( db, "by_site", by_site, append = TRUE, overwrite = FALSE )
	DBI::dbExecute( db, "DELETE FROM by_sample WHERE source == 'MalariaGEN Pf7' ")
	DBI::dbWriteTable( db, "by_sample", by_sample, append = TRUE, overwrite = FALSE )
	DBI::dbDisconnect( db )
