library( RSQLite )
library( dplyr )
db = dbConnect( dbDriver( "SQLite" ), "input/hbs-pf-pf8.sqlite" )
D = dbGetQuery( db, "SELECT * FROM by_sample" )

# Counts for Supplementary Table 2
loci = c( "Pfsa1", "Pfsa2", "Pfsa3", "Pfsa4")

summary = (
	D
	%>% filter( exclude == 'no' & !is.na( ref ) & locus %in% loci )
	%>% group_by( source, datatype, country, year, site, latitude, longitude, ID )
	%>% summarise( n_records = n() )
	%>% group_by( source, country, datatype )
	%>% summarise(
		mean_latitude = mean( latitude ),
		mean_longitude = mean( longitude ),
		number_of_samples = n(),
		number_of_sites = length( unique( sprintf( "%.6f-%.6f", longitude, latitude )))
	)
)

locus_summaries = list()
for( l in loci ) {
	locus_summaries[[l]] = (
		D
		%>% filter( exclude == 'no' & !is.na( ref ) & locus == l )
		%>% group_by( source, datatype, country, year, site, latitude, longitude, ID, ref, mixed, nonref )
		%>% summarise( n_records = n() )
		%>% group_by( source, country, datatype )
		%>% summarise(
			ref = sum( ref ),
			mixed = sum( mixed ),
			nonref = sum( nonref )
		)
	)
	colnames(locus_summaries[[l]])[4:6] = sprintf( "%s:%s", l, colnames(locus_summaries[[l]])[4:6] )
}
sum( summary$number_of_samples )
dim(summary)
summary = (
	summary
	%>% left_join( locus_summaries$`Pfsa1`, by = c( "source", "country", "datatype" ))
	%>% left_join( locus_summaries$`Pfsa2`, by = c( "source", "country", "datatype" ))
	%>% left_join( locus_summaries$`Pfsa3`, by = c( "source", "country", "datatype" ))
	%>% left_join( locus_summaries$`Pfsa4`, by = c( "source", "country", "datatype" ))
	%>% arrange( mean_longitude )
)
dim(summary)

dir.create( "output/pf=pf8-version/tables", showWarnings = FALSE )
readr::write_csv(
	summary
	%>% arrange( mean_longitude )
	%>% mutate( blank1 = '', blank2 = '', blank3 = '', blank4 = '' )
	%>% select(
		Source = source,
		Country = country,
		`Data type` = datatype,
		`Number of samples` = number_of_samples,
		`Number of sites` = number_of_sites,
		blank1,
		ref1   = `Pfsa1:ref`,
		mixed1 = `Pfsa1:mixed`,
		alt1   = `Pfsa1:nonref`,
		blank2,
		ref2   = `Pfsa2:ref`,
		mixed2 = `Pfsa2:mixed`,
		alt2   = `Pfsa2:nonref`,
		blank3,
		ref3   = `Pfsa3:ref`,
		mixed3 = `Pfsa3:mixed`,
		alt3   = `Pfsa3:nonref`,
		blank4,
		ref4   = `Pfsa4:ref`,
		mixed4 = `Pfsa4:mixed`,
		alt4   = `Pfsa4:nonref`
	),
	file = "output/pf=pf8-version/tables/table_S2_pf_counts.csv"
)

# Pfsa1 counts, to match table in ,figure 1
summary = (
	D
	%>% filter( exclude == 'no' & !is.na( ref ) & (mixed == 0) & locus %in% c( "Pfsa1" ))
	%>% group_by( source, datatype, country, year, site, latitude, longitude, ID )
	%>% summarise( n_records = n() )
	%>% group_by( source, country, datatype )
	%>% summarise(
		mean_latitude = mean( latitude ),
		mean_longitude = mean( longitude ),
		number_of_samples = n(),
		number_of_sites = length( unique( sprintf( "%.4f-%.4f", longitude, latitude )))
	)
)
sum( summary$number_of_samples )

readr::write_csv( summary, file = "output/pf=pf8-version/tables/figure1_pf_counts.csv" )


# site count to match Fig 1 legend
summary = (
	D
	%>% filter( exclude == 'no' & !is.na( ref ) & locus %in% c( "Pfsa1", "Pfsa2", "Pfsa3", "Pfsa4" ))
	%>% group_by( country, latitude, longitude )
	%>% summarise( number_of_samples = n() )
	%>% group_by( country )
	%>% summarise( number_of_sites = n() )
)

# Check Pfsa2 / Pfsa4 sample counts
db = dbConnect( dbDriver( "SQLite" ), "input/hbs-pf-pf8.sqlite" )
D = dbGetQuery( db, "SELECT * FROM by_sample" )
summary = (
	D
	%>% filter( exclude == 'no' & !is.na( ref ) & (mixed==0) & locus %in% c( "Pfsa2" ))
	%>% group_by( source, datatype, country, year, site, latitude, longitude, ID, ref, mixed, nonref )
	%>% summarise( n_records = n() )
	%>% group_by( country )
	%>% summarise(
		ref = sum( ref ),
		mixed = sum( mixed ),
		nonref = sum( nonref ),
		mean_latitude = mean( latitude ),
		mean_longitude = mean( longitude ),
		number_of_samples = n(),
		number_of_sites = length( unique( sprintf( "%.6f-%.6f", longitude, latitude )))
	)
	%>% mutate(
		n = ref + nonref,
		`f+` = nonref/n
	)
)
print( summary %>% filter( `f+` > 0.01), n = 100, width = 1000 )


summary = (
	D
	%>% filter( exclude == 'no' & !is.na( ref ) & (mixed==0) & locus %in% c( "Pfsa4" ))
	%>% group_by( source, datatype, country, year, site, latitude, longitude, ID, ref, mixed, nonref )
	%>% summarise( n_records = n() )
	%>% group_by( country )
	%>% summarise(
		ref = sum( ref ),
		mixed = sum( mixed ),
		nonref = sum( nonref ),
		mean_latitude = mean( latitude ),
		mean_longitude = mean( longitude ),
		number_of_samples = n(),
		number_of_sites = length( unique( sprintf( "%.6f-%.6f", longitude, latitude )))
	)
	%>% mutate(
		n = ref + nonref,
		`f+` = ref/n
	)
)
print( summary %>% filter( `f+` > 0.01), n = 100, width = 1000 )
