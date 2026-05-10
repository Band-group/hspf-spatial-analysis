library( tidyverse )
library( rbgen )
library( argparse )

source("scripts/load.plasmodb.genes.R")

echo <- function( message, ... ) {
	cat( sprintf( message, ... ))
}

blank.plot <- function( xlim = c(0,1), ylim = c(0,1), ... ) {
	plot( 0, 0, col = 'white', bty = 'n', xaxt = 'n', yaxt = 'n', xlim = xlim, ylim = ylim, ... )
}

arg_sets = list(
	`Pfsa1` = list(
		pf7 = "outputs/pf7/vcf/07_ancestral/Pf3D7_02_v3.bgen",
		samples = "outputs/pf7/samples/filtered_samples.tsv",
		genes = "genes/PlasmoDB-65_Pfalciparum3D7.gff.gz",
		gaf = "genes/PlasmoDB-65_Pfalciparum3D7_GO.gaf.gz",
		trees = c(
			`Pfsa1` = "outputs/pf7/relate/output/trees/popsize/pf7.relate.Pf3D7_02_v3.Ne=100000.mu=6.3162e-9.dpg=365.bp=631190.newick",
			`Pfsa2` = "outputs/pf7/relate/output/trees/popsize/pf7.relate.Pf3D7_02_v3.Ne=100000.mu=6.3162e-9.dpg=365.bp=814288.newick",
			`Pfsa3` = "outputs/pf7/relate/output/trees/popsize/pf7.relate.Pf3D7_11_v3.Ne=100000.mu=6.3162e-9.dpg=365.bp=1058035.newick",
			`Pfsa4` = "outputs/pf7/relate/output/trees/popsize/pf7.relate.Pf3D7_04_v3.Ne=100000.mu=6.3162e-9.dpg=365.bp=1121472.newick",
			`CRT`   = "outputs/pf7/relate/output/trees/popsize/pf7.relate.Pf3D7_07_v3.Ne=100000.mu=6.3162e-9.dpg=365.bp=403625.newick"
		),
		length_samples = c(
			`Pfsa1` = "outputs/pf7/relate/output/trees/popsize/pf7.relate.Pf3D7_02_v3.Ne=100000.mu=6.3162e-9.dpg=365.bp=631190.samples.newick"
		),
		annotation = "outputs/pf7/vcf/04_merged/Pf3D7_02_v3.merged.annotation.tsv.gz",
		margin = 20000,
		min_maf = 0.005,
		focus_margin = 10000,
		focus = "Pf3D7_02_v3:631190",
		zoom_region = list(
			chromosome = "Pf3D7_02_v3",
			start = 626250,
			end   = 633750
		),
		countries = NULL,
		tree_annotated_positions = 631190
	),
	`Pfsa3` = list(
		pf7 = "outputs/pf7/vcf/07_ancestral/Pf3D7_11_v3.bgen",
		samples = "outputs/pf7/samples/filtered_samples.tsv",
		genes = "genes/PlasmoDB-65_Pfalciparum3D7.gff.gz",
		gaf = "genes/PlasmoDB-65_Pfalciparum3D7_GO.gaf.gz",
		trees = c(
			`Pfsa3` = "outputs/pf7/relate/output/trees/popsize/pf7.relate.Pf3D7_11_v3.Ne=100000.mu=6.3162e-9.dpg=365.bp=1058035.newick"
		),
		length_samples = c(
			`Pfsa3` = "outputs/pf7/relate/output/trees/popsize/pf7.relate.Pf3D7_11_v3.Ne=100000.mu=6.3162e-9.dpg=365.bp=1058035.samples.newick"
		),
		annotation = "outputs/pf7/vcf/04_merged/Pf3D7_11_v3.merged.annotation.tsv.gz",
		margin = 20000,
		min_maf = 0.005,
		focus_margin = 10000,
		focus = "Pf3D7_11_v3:1057437",
		zoom_region = list(
			chromosome = "Pf3D7_11_v3",
			start = 1053000,
			end   = 1060000
		),
		countries = NULL
	),
	`Pfsa2` = list(
		pf7 = "outputs/pf7/vcf/07_ancestral/Pf3D7_02_v3.bgen",
		samples = "outputs/pf7/samples/filtered_samples.tsv",
		genes = "genes/PlasmoDB-65_Pfalciparum3D7.gff.gz",
		gaf = "genes/PlasmoDB-65_Pfalciparum3D7_GO.gaf.gz",
		trees = c(
			`Pfsa2` = "outputs/pf7/relate/output/trees/popsize/pf7.relate.Pf3D7_02_v3.Ne=100000.mu=6.3162e-9.dpg=365.bp=814288.newick"
		),
		length_samples = c(
			`Pfsa2` = "outputs/pf7/relate/output/trees/popsize/pf7.relate.Pf3D7_02_v3.Ne=100000.mu=6.3162e-9.dpg=365.bp=814288.samples.newick"
		),
		annotation = "outputs/pf7/vcf/04_merged/Pf3D7_02_v3.merged.annotation.tsv.gz",
		margin = 20000,
		min_maf = 0.005,
		focus_margin = 10000,
		focus = "Pf3D7_02_v3:814288",
		zoom_region = list(
			chromosome = "Pf3D7_02_v3",
			start = 811288,
			end   = 817288
		),
		countries = NULL,
		tree_annotate_all = TRUE
#			'Democratic_Republic_of_the_Congo', 'Kenya', 'Rwanda', 'Uganda', 'Malawi', 'Zambia', 'Mozambique', 'Tanzania'
	),
	`Pfsa4` = list(
		pf7 = "outputs/pf7/vcf/07_ancestral/Pf3D7_04_v3.bgen",
		samples = "outputs/pf7/samples/filtered_samples.tsv",
		genes = "genes/PlasmoDB-65_Pfalciparum3D7.gff.gz",
		gaf = "genes/PlasmoDB-65_Pfalciparum3D7_GO.gaf.gz",
		trees = c(
			`Pfsa4` = "outputs/pf7/relate/output/trees/popsize/pf7.relate.Pf3D7_04_v3.Ne=100000.mu=6.3162e-9.dpg=365.bp=1121472.newick"
		),
		length_samples = c(
			`Pfsa4` = "outputs/pf7/relate/output/trees/popsize/pf7.relate.Pf3D7_04_v3.Ne=100000.mu=6.3162e-9.dpg=365.bp=1121472.samples.newick"
		),
		annotation = "outputs/pf7/vcf/04_merged/Pf3D7_04_v3.merged.annotation.tsv.gz",
		margin = 20000,
		min_maf = 0.005,
		focus_margin = 10000,
		focus = "Pf3D7_04_v3:1121472",
		zoom_region = list(
			chromosome = "Pf3D7_04_v3",
			start = 1116472,
			end   = 1126472
		),
		countries = c('Mauritania', 'Senegal', 'Gambia', 'Guinea', 'Mali', 'Burkina_Faso', 'Cote_dIvoire', 'Ghana', 'Togo', 'Benin', 'Nigeria', 'Cameroon', 'Gabon'),
		tree_annotate_all = FALSE
	)
)

#locus = "Pfsa1"
#locus = "Pfsa2"
#locus = "Pfsa3"
locus = "Pfsa4"
args = arg_sets[[locus]]

focus = tibble(
	chromosome = strsplit( args$focus, split = ':' )[[1]][1],
	position = as.integer( strsplit( args$focus, split = ':' )[[1]][2] )
)
focus$start = focus$position - args$margin
focus$end = focus$position + args$margin

#args = parse_arguments()
stopifnot( args$focus_margin <= args$margin )
echo( "++ Loading samples from %s...\n", args$samples )

regions = c(
	Gambia = "west",
	Senegal = "west",
	Guinea = "west",
	Mauritania = "west",
	Cote_dIvoire = "west",
	Mali = "west",
	Burkina_Faso = "west",
	Ghana = "west",
	Benin = "west",
	Nigeria = "west",
	Gabon = "west",
	Cameroon = "west",
	Democratic_Republic_of_the_Congo = "central",
	Sudan = "east",
	Uganda = "east",
	Malawi = "east",
	Tanzania = "east",
	Mozambique = "east",
	Kenya = "east",
	Ethiopia = "east",
	Madagascar = "east"
)
samples = (
	read_tsv( args$samples )
	%>% mutate( relate_sample_index = sprintf( "%d", 0:(length(Sample)-1)))
	%>% mutate( region = regions[Country] )
)
echo( "++ Ok, %d samples loaded.\n", nrow( samples ))

echo( "++ focussing on region:\n" )
print(focus)

# Turn pop and country into factors
populations = unique( samples$Population[ order( samples$`Country longitude`)] )
countries = unique( (samples %>% arrange( `Country longitude`))$Country )
samples$Population = factor( samples$Population, levels = populations )
samples$Country = factor( samples$Country, levels = countries )
samples$region = factor( samples$region, levels = c( "west", "central", "east") )

if( is.null(args$countries)) {
	args$countries = levels( samples$Country )
}

H = load.genotypes( args$pf7, focus )

# The data is haploid, and usually biallelic, just take the alt / 2nd allele calls 
HD = H$data[,,2]
variants = as_tibble(H$variants)
variants$freq = rowSums( HD, na.rm = T ) / rowSums( !is.na( HD ))
annotations = simplify_alleles(
	variants %>% select( chromosome, position )
	%>% inner_join(
		readr::read_tsv( args$annotation ),
		by = c( "chromosome", "position" )
	)
)
variants = bind_cols(
	variants,
	purrr::map_dfr( 1:nrow(variants), function(i) { split_annotations( variants[i,], annotations$annotation[i] ) })
)
variants$consequence[ variants$consequence_allele != variants$allele0 & variants$consequence_allele != variants$allele1 ] = 'none'

focus.variant = which( variants$position == focus$position )
wIn = which(
	variants$freq >= args$min_maf
	& variants$freq <= (1-args$min_maf)
	& variants$position >= focus$position - 10000
	& variants$position <= focus$position + 10000
)

samples$focus_genotype = HD[ which( variants$position == focus$position), ]

{
	spec = list(
		locus = locus,
		focus = focus,
		zoom_region = args$zoom_region,
		countries = args$countries,
		variants = variants[wIn,],
		samples = (
			samples
			%>% mutate( index = sprintf( "%d", 1:nrow( samples )))
			%>% filter( Country %in% args$countries )
#			%>% group_by( Country )
			%>% group_by( Country, focus_genotype )
			%>% slice_sample( n = 20 )
		),
		genes = (
			load.plasmodb.genes( gff = args$genes, gaf = args$gaf )
			%>% filter( start <= focus$end & end >= focus$start )
			%>% filter(
				(ID == 'PF3D7_0215300' | Parent == 'PF3D7_0215300' | Parent == 'PF3D7_0215300.1')
				| (Parent == 'PF3D7_1127000' | ID == 'PF3D7_1127000' | Parent == 'PF3D7_1127000.1')
				| (Parent == 'PF3D7_1126900' | ID == 'PF3D7_1126900' | Parent == 'PF3D7_1126900.1')
				| (ID == 'PF3D7_0220300' | Parent == 'PF3D7_0220300' | Parent == 'PF3D7_0220300.1')
				| (ID == 'PF3D7_0424700' | Parent == 'PF3D7_0424700' | Parent == 'PF3D7_0424700.1')
			)
		)
	)
	spec$haplotypes = HD[wIn, spec$samples$Sample]
	spec$annotated_variants = (
		(
			find_high_ld_variants( HD, which( variants$position == focus$position ))
			%>% filter( position >= focus$position - 5000 & position <= focus$position + 5000 )
			%>% mutate( `f+` = (`++`/(`++`+`+-`)), `f-` = (`-+`/(`-+`+`--`)) )
			%>% filter( `f-` < 0.05 & (`f+`/`f-`) > 5 & freq > 0.02 )
			%>% filter( position %in% spec$variants$position )
			#%>% filter( freq >= 0.02 )
			%>% mutate(
				shape = 21,
				size = 0.75,
				text.size = 0.5,
				colour = consequence.colours[ consequence ],
				border = 'black',
				font = 1
			)
		)
	)
	spec$annotated_variants$shape[ spec$annotated_variants$position == focus$position ] = 25
	spec$annotated_variants$size[ spec$annotated_variants$position == focus$position ] = 1.25
	spec$annotated_variants$text.size[ spec$annotated_variants$position == focus$position ] = 1
	spec$annotated_variants$font[ spec$annotated_variants$position == focus$position ] = 1

	# Load tree
	spec$trees = list()
	spec$fulltrees = list()
	for( name in names( args$trees )) {
		tree = ape::read.tree( args$trees[[name]] )
		# This bit gets the sub-tree for the chosen samples.
		spec$fulltrees[[name]] = tree
		spec$fulltrees[[name]]$tip.sample = samples$Sample[ match( spec$fulltrees[[name]]$tip.label, samples$relate_sample_index )]
		spec$trees[[name]] = ape::keep.tip( tree, spec$samples$relate_sample_index )
		spec$trees[[name]]$tip.sample = samples$Sample[ match( spec$trees[[name]]$tip.label, samples$relate_sample_index )]
	}
	if( args$tree_annotate_all ) {
		spec$tree_mutations = assign.mutations(
			spec$trees[[locus]],
			spec$annotated_variants,
			spec$haplotypes[ spec$variants$position %in% spec$annotated_variants$position,, drop = F ],
			threshold = 0.99
		)
	} else {
		spec$tree_mutations = assign.mutations(
			spec$trees[[locus]],
			spec$annotated_variants %>% filter( position == focus$position ),
			spec$haplotypes[ spec$variants$position == focus$position,, drop = F ],
			threshold = 0.99
		)
	}

	spec$fulltree_mutations = assign.mutations(
		spec$fulltrees[[locus]],
		spec$annotated_variants %>% filter( position == focus$position ),
		HD[variants$position == focus$position,, drop = F ],
		verbose = TRUE
	)

	spec$length_samples = list()
	for( name in names( args$length_samples )) {
		spec$length_samples[[name]] = list()
		X = readr::read_tsv( args$length_samples[[name]] )
		spec$length_samples[[name]] = lapply(
			1:nrow(X),
			function(i) {
				ape::read.tree( text = X$tree[i] )
			}
		)
	}
	# Compute age range estimates from samples
	{
		m = spec$fulltree_mutations[1,]
		upperlower = sapply(
			1:length( spec$length_samples[[locus]] ),
			function(i) {
				tree = spec$length_samples[[locus]][[i]]
				times = ape::node.depth.edgelength( tree )
				tmcra = max( times )
				time.ago = tmcra - times
				return( c(
					lower = time.ago[m$node],
					upper = time.ago[m$parent]
				)) ;
			}
		)
		spec$age_range = tibble::tibble(
			lower2.5 = quantile( upperlower['lower',], 0.025 ),
			upper97.5 = quantile( upperlower['upper',], 0.975 ),
			lower25 = quantile( upperlower['lower',], 0.25 ),
			upper75 = quantile( upperlower['upper',], 0.75 )
		)
	}


	# Fix sample ordering for tree
	ho = match( spec$trees[[locus]]$tip.sample, spec$samples$Sample )
	stopifnot( length( which( is.na( ho ))) == 0 )
	spec$samples = spec$samples[ho,]
	spec$haplotypes = spec$haplotypes[,ho]
}

{
	source( "../spatial/code/functions.R" )
	source( "scripts/layout.intervals.R" )
	source( "scripts/plot.genes.R" )
	source( "scripts/haplotype_figure_impl.R" )

	figure_3(
		spec = spec,
		colour.column = "Country",
		split = c( 0.425, 0.575 ),
		width = 10,
		height = 6,
		sprintf( "outputs/figures/figure_3-%s:%d.pdf", focus$chromosome, focus$position ),
		colours = list(
			Country = country.colours()[spec$countries],
			region = c(
				west = country.colours()[["Gambia"]],
				central = country.colours()[["Democratic_Republic_of_the_Congo"]],
				east = country.colours()[["Kenya"]]
			)
		)
	)	
}
