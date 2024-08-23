library( tidyverse )
library( dplyr )
library( rbgen )

samples = read_tsv( '/well/band/projects/pf7/data/samples/Pf7_samples.txt' )
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

G = bgen.load(
	"/well/band/projects/pf7/results/bgen/pf7.filtered.bgen",
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

for( i in 1:nrow( G$variants )) {
	samples[,G$variants$ID[i]] = G$dosage[i,]
}

wIn = which( samples$`Exclusion reason` == 'Analysis_set' )
write_tsv( samples, "results/pfsa_genotypes/pf7.tsv" )

by_site_and_year = (
	samples
	%>% filter(
		`Exclusion reason` == 'Analysis_set'
	)
	%>% group_by(
		Study, Country, `Admin level 1`, `Admin level 1 latitude`, `Admin level 1 longitude`, `Year`
	) %>% summarise(
		`Pfsa1:ref-` = sum( !is.na( `chr2:631190:T>A` )) - sum( `chr2:631190:T>A`, na.rm = T  ),
		`Pfsa1:nonref` = sum( `chr2:631190:T>A`, na.rm = T  ),
		`Pfsa2:ref` = sum( !is.na( `chr2:814288:C>T` )) - sum( `chr2:814288:C>T`, na.rm = T  ),
		`Pfsa2:nonref` = sum( `chr2:814288:C>T`, na.rm = T  ),
		`Pfsa3:ref` = sum( !is.na( `chr11:1058035:T>A` )) - sum( `chr11:1058035:T>A`, na.rm = T  ),
		`Pfsa3:nonref` = sum( `chr11:1058035:T>A`, na.rm = T  ),
		`Pfsa4:ref` = sum( !is.na( `chr4:1121472:T>A` )) - sum( `chr4:1121472:T>A`, na.rm = T  ),
		`Pfsa4:nonref` = sum( `chr4:1121472:T>A`, na.rm = T  )
	)
)

write_tsv(
	by_site_and_year,
	"resuilts/genotypes/pf7_by_site_and_year.tsv"
)

by_site = (
	samples
	%>% filter(
		`Exclusion reason` == 'Analysis_set'
	)
	%>% group_by(
		Study, Country, `Admin level 1`, `Admin level 1 latitude`, `Admin level 1 longitude`
	) %>% summarise(
		N = n(),
		`Pfsa1:ref` = sum( !is.na( `chr2:631190:T>A` )) - sum( `chr2:631190:T>A`, na.rm = T ),
		`Pfsa1:nonref` = sum( `chr2:631190:T>A`, na.rm = T ),
		`Pfsa2:ref` = sum( !is.na( `chr2:814288:C>T` )) - sum( `chr2:814288:C>T`, na.rm = T ),
		`Pfsa2:nonref` = sum( `chr2:814288:C>T`, na.rm = T ),
		`Pfsa3:ref` = sum( !is.na( `chr11:1058035:T>A` )) - sum( `chr11:1058035:T>A`, na.rm = T ),
		`Pfsa3:nonref` = sum( `chr11:1058035:T>A`, na.rm = T ),
		`Pfsa4:ref` = sum( !is.na( `chr4:1121472:T>A` )) - sum( `chr4:1121472:T>A`, na.rm = T ),
		`Pfsa4:nonref` = sum( `chr4:1121472:T>A`, na.rm = T )
	)
)
by_site = bind_cols(
	source = "MalariaGEN Pf7",
	by_site
)
colnames( by_site )[1:6] = c(
	"source", "study", "country", "site", "latitude", "longitude"
)


write_tsv(
	by_site,
	"results/genotypes/pf7_by_site.tsv.gz"
)

db = DBI::dbConnect( RSQLite::SQLite(), "results/genotypes/hbs-pf.sqlite" )
DBI::dbWriteTable( db, "by_site", by_site, overwrite = TRUE )
DBI::dbDisconnect( db )
