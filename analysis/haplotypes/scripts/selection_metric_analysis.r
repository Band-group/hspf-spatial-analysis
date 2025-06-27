ihs = readr::read_tsv( "outputs/pf7/selscan/output/pf7.selscan.ihs.bins=2.5%.tsv.gz" )
ihh12 = readr::read_tsv( "outputs/pf7/selscan/output/pf7.selscan.ihh12.bins=2.5%.tsv.gz" )
beta = readr::read_tsv( "outputs/pf7/betascan/advanced/pf7.betascan.window=5000.p=50.tsv.gz" )
source( 'scripts/load.plasmodb.genes.R' )
source( 'scripts/annotate.nearest.genes.R' )
gff = load.plasmodb.genes(
	"../../../../pfsa/data/genes/pf/3D7/PlasmoDB-65_Pfalciparum3D7.gff.gz",
	"../../../../pfsa/data/genes/pf/3D7/PlasmoDB-65_Pfalciparum3D7_GO.gaf.gz"
)
genes = (
	gff
	%>% filter( type %in% c( 'protein_coding_gene', 'pseudogene' ))
)
genes$known_antigenic = 0
genes$known_antigenic[
	grep(
		"rifin|stevor|PfEMP|merozoite surface|erythrocyte surface|surface-associated|erythrocyte binding|reticulocyte binding|antigen|cytoadherence",
		genes$attributes
	)
] = 1
positions = unique(
	rbind(
		ihs[,c("chromosome", "position")],
		ihh12[,c("chromosome", "position")],
		beta[,c("chromosome", "position")]
	)
)
positions = bind_cols(
	positions,
	annotate.nearest.genes( positions, genes, margin = 5000 )
)
positions$region_known_antigenic = sapply(
	stringr::str_split( positions$genes_in_region, ";" ),
	function(s) { 
		M = match( s, genes$ID )
		return( max( genes$known_antigenic[M] ))
	}
)

combined = (
	positions
	%>% full_join(
		ihs %>% filter( years == 'all' ) %>% select( country, chromosome, position, frequency, uIHS ),
		by = c( "chromosome", "position" )
	)
	%>% full_join(
		ihh12 %>% filter( years == 'all' ) %>% select( country, chromosome, position, uiHH12 ),
		by = c( "country", "chromosome", "position" )
	)
	%>% full_join(
		beta %>% select(
			country, chromosome, position, derived, total, variants_in_window,
			beta, tajimas_d, dango, pi,
			pi_fraction_between
		),
		by = c( "country", "chromosome", "position" )
	)
	%>% filter( region_known_antigenic == 0 )
	%>% tidyr::pivot_longer(
		cols = c( "uIHS", "uiHH12", "beta", "tajimas_d", "dango", "pi", "pi_fraction_between" ),
		names_to = "statistic"
	)
)

focus = combined %>% filter( position %in% c( 631190, 814288, 1058035, 1057437, 1121472 ))
m = 0.01
comparison = (
	focus %>% select(
		country,
		focus_chromosome = chromosome, focus_position = position, focus_frequency = frequency,
		derived, total, variants_in_window, statistic,
		focus_value = value
	)
	%>% inner_join(
		combined %>% select( country, frequency, statistic, value ),
		by = c( "country", "statistic" ),
		relationship = "many-to-many"
	)
	%>% filter(
		(frequency >= focus_frequency - m) & (frequency <= focus_frequency + m)
	)
	%>% group_by( country, focus_chromosome, focus_position, focus_frequency, derived, total, statistic, focus_value )
	%>% summarise(
		n = sum( !is.na( value )),
		n_above = sum( value >= focus_value ),
		`value:mean` = mean( value, na.rm = T ),
		`value:sd` = sd( value, na.rm = T )
	)
	%>% mutate(
		`value:norm` = (focus_value - `value:mean`) / `value:sd`
	)
	%>% select(
		country,
		chromosome = focus_chromosome, position = focus_position, frequency = focus_frequency,
		derived, total,
		statistic,
		value = focus_value,
		comparison_n = n,
		comparison_n_above = n_above,
		`value:mean`,
		`value:sd`,
		`value:norm`
	)
	%>% mutate(
		comparison_propn_above = comparison_n_above / comparison_n
	)
)

readr::write_csv( comparison, file = "outputs/pf7/selection.csv" )
