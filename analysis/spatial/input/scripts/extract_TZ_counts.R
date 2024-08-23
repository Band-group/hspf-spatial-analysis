library( tidyverse )
library( dplyr )
library( dbplyr )
library( rbgen )

samples = read_tsv( 'data/tanzania/Moser_et_al_2021/all_metadata_header.txt' )
controls = scan( 'data/tanzania/Moser_et_al_2021/list_of_controls', what = character() ) 
samples$is_control = 0
samples$is_control[ samples$ID %in% controls ] = 1
# These lat/long values are taken from Google Maps
# in comparison to the map in Figure 1 of Moser et al 2021 https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8088766/
latlong = list(
	"UVINZA" = c( -5.105266, 30.383864 ),
	"BUHIGWE" = c( -4.497555, 29.894095 ),
	"CHATO" = c( -2.638419, 31.768127 ),
	"NYANGHWALE" = c( -3.196338, 32.649443 ),
	"NYASA" = c( -11.378040, 35.126810 ),
	"TUNDURU" = c( -11.040612, 37.333045 ),
	"NANYUMBU" = c( -11.144070, 38.492533 ),
	"MTWARA DC" = c( -10.307817, 40.177817 ),
	"KIBAHA" = c( -6.781600, 38.990543 ),
	"KIGOMA-UJIJI" = c( -4.911249, 29.674744 ),
	"ILEMELA" = c( -2.515741, 32.909151 ),
	"KYELA" = c( -9.590593, 33.867588 ),
	"MASASI" = c( -10.729617, 38.806608 )
)
samples$latitude = sapply(
	samples$DISTRICT, function(d) { latlong[[d]][1] ; }
)
samples$longitude = sapply(
	samples$DISTRICT, function(d) { latlong[[d]][2] ; }
)

SNPs = tibble(
	chromosome = sprintf( "chr%d", c( 2, 2, 2, 2, 2, 4, 4, 11, 11 )),
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
	"data/tanzania/Moser_et_al_2021/IBC_variants.fixed_genos.biallelic.targets_only.recode.bgen",
	ranges = tibble(
		chromosome = SNPs$chromosome,
		start = SNPs$position,
		end = SNPs$position
	)
)
# for the sake of this analysis, we ignore mixed calls
# or calls of other allels
G$data = G$data[,,c(1,3)]
G$dosage = G$data[,,2]
G$dosage[ rowSums(G$data,dims=2) == 0 ] = NA

G$variants$ID = sprintf(
	"chr%s:%d:%s>%s",
	gsub( "chr", "", G$variants$chromosome ),
	G$variants$position,
	G$variants$allele0,
	G$variants$allele1
)

samples = samples[ match( G$samples, samples$sample_id ), ]
for( i in 1:nrow( G$variants )) {
	samples[,G$variants$ID[i]] = G$dosage[i,]
}

write_tsv( samples, "results/pfsa_genotypes/Moser_et_al_2021.tsv" )

by_site = (
	samples
	%>% filter( is_control == 0 )
	%>% group_by(
		DISTRICT, REGION, latitude, longitude
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
by_site$source = "Moser et al 2021"
by_site$study = "Moser et al 2021"
by_site$country = "Tanzania"
by_site$site = sprintf( "%s:%s", by_site$REGION, by_site$DISTRICT )
by_site = by_site[,
	c(
		"source", "study", "country", "site", "latitude", "longitude",
		"N",
		"Pfsa1:ref", "Pfsa1:nonref",
		"Pfsa2:ref", "Pfsa2:nonref",
		"Pfsa3:ref", "Pfsa3:nonref",
		"Pfsa4:ref", "Pfsa4:nonref"
	)
]

write_tsv( samples, "results/genotypes/Moser_et_al_2021_by_site.tsv" )

db = DBI::dbConnect( RSQLite::SQLite(), "results/genotypes/hbs-pf.sqlite" )
DBI::dbExecute( db, "DELETE FROM by_site WHERE source == 'Moser_et_al_2021' ")
DBI::dbWriteTable( db, "by_site", by_site, overwrite = FALSE, append = TRUE )
DBI::dbDisconnect( db )
