library( dplyr )
library( argparse )

options(width=300)
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
		default = "../../../data/uganda"
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
	data = sprintf( "%s/pfsa_data_uganda_wgs.tsv", args$indir )
)

X = (
	readr::read_tsv( paths$data )
)

X$Pf3D7_02_v3_631190_T_A = c(
	'T|T' = 0,
	'T/T' = 0,
	'T/A' = NA,
	'T|A' = NA,
	'A/A' = 1,
	'A|A' = 1
)[ X$Pf3D7_02_v3_631190_T_A ]

X$Pf3D7_02_v3_814288_C_T = c(
	'C/C' = 0,
	'C|C' = 0,
	'C/T' = NA,
	'C|T' = NA,
	'T/T' = 1,
	'T|T' = 1
)[ X$Pf3D7_02_v3_814288_C_T ]

X$Pf3D7_04_v3_1121472_T_A = c(
	'T|T' = 0,
	'T/T' = 0,
	'T/A' = NA,
	'T|A' = NA,
	'A/A' = 1,
	'A|A' = 1
)[ X$Pf3D7_04_v3_1121472_T_A ]

X$Pf3D7_11_v3_1058035_T_A = c(
	'T|T' = 0,
	'T/T' = 0,
	'T/A' = NA,
	'T|A' = NA,
	'A/A' = 1,
	'A|A' = 1
)[ X$Pf3D7_11_v3_1058035_T_A ]

echo( "Creating by_sample...\n" )
by_sample = (
	X
	%>% mutate(
		source = "Greenwood Uganda 2017-2022",
		study = "Greenwood Uganda 2017-2022",
		datatype = 'WGS',
		country = "Uganda",
		site,
		latitude,
		longitude,
		ID = sample_name, 
		N = 1,
		`Pfsa1:ref`    = 1 - Pf3D7_02_v3_631190_T_A,
		`Pfsa1:nonref` = Pf3D7_02_v3_631190_T_A,
		`Pfsa2:ref`    = 1 - Pf3D7_02_v3_814288_C_T,
		`Pfsa2:nonref` = Pf3D7_02_v3_814288_C_T,
		`Pfsa3:ref` = 1 - Pf3D7_11_v3_1058035_T_A,
		`Pfsa3:nonref` = Pf3D7_11_v3_1058035_T_A,
		`Pfsa4:ref` = 1 - Pf3D7_04_v3_1121472_T_A,
		`Pfsa4:nonref` = Pf3D7_04_v3_1121472_T_A,
		exclude = "no"
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
print( by_sample )
echo( "Creating by_site...\n" )

by_site = (
	by_sample
	%>% filter( !is.na( latitude ))
	%>% filter( exclude == 'no' )
	%>% group_by(
		source, study, datatype, country, site, latitude, longitude
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

echo( "Writing to db...\n" )

db = DBI::dbConnect( RSQLite::SQLite(), args$output )
DBI::dbExecute( db, "DELETE FROM by_site WHERE source == 'Greenwood_Uganda_2017-2022' ")
DBI::dbWriteTable( db, "by_site", by_site, overwrite = FALSE, append = TRUE )
DBI::dbExecute( db, "DELETE FROM by_sample WHERE source == 'Greenwood_Uganda_2017-2022' ")
DBI::dbWriteTable( db, "by_sample", by_sample, overwrite = FALSE, append = TRUE )
DBI::dbDisconnect( db )
