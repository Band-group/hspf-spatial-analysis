library( tidyverse )
library( rbgen )
library( argparse )

source("scripts/load.plasmodb.genes.R")
source( "scripts/haplotype_figure_impl.R" )
source( "scripts/figure4_impl.R" )

echo <- function( message, ... ) {
	cat( sprintf( message, ... ))
}

blank.plot <- function( xlim = c(0,1), ylim = c(0,1), ... ) {
	plot( 0, 0, col = 'white', bty = 'n', xaxt = 'n', yaxt = 'n', xlim = xlim, ylim = ylim, ... )
}

args = list(
	min_maf = 0.005,
	countries = NULL,
	samples   = "outputs/pf7/samples/filtered_samples.tsv",
	genes    = "genes/PlasmoDB-65_Pfalciparum3D7.gff.gz",
	gaf      = "genes/PlasmoDB-65_Pfalciparum3D7_GO.gaf.gz",
	regions = list(
		`Pfsa1` = list(
			tree        = "outputs/pf7/relate/output/trees/popsize/pf7.relate.Pf3D7_02_v3.Ne=100000.mu=6.3162e-9.dpg=365.bp=631190.newick",
			treesamples = "outputs/pf7/relate/output/trees/popsize/pf7.relate.Pf3D7_02_v3.Ne=100000.mu=6.3162e-9.dpg=365.bp=631190.samples.newick",
			genotypes   = "outputs/pf7/vcf/07_ancestral/Pf3D7_02_v3.bgen",
			annotation  = "outputs/pf7/vcf/04_merged/Pf3D7_02_v3.merged.annotation.tsv.gz",
			focus       = tibble::tibble(
				chromosome  = 'Pf3D7_02_v3',
				position    = 631190,
				start       = 621190,
				end         = 641190,
				zoom_start  = 626250,
				zoom_end    = 633750
			),
			tree_annotated_positions = 631190
		),
		`Pfsa3` = list(
			tree        = "outputs/pf7/relate/output/trees/popsize/pf7.relate.Pf3D7_11_v3.Ne=100000.mu=6.3162e-9.dpg=365.bp=1058035.newick",
			treesamples = "outputs/pf7/relate/output/trees/popsize/pf7.relate.Pf3D7_11_v3.Ne=100000.mu=6.3162e-9.dpg=365.bp=1058035.samples.newick",
			genotypes   = "outputs/pf7/vcf/07_ancestral/Pf3D7_11_v3.bgen",
			annotation  = "outputs/pf7/vcf/04_merged/Pf3D7_11_v3.merged.annotation.tsv.gz",
			focus       = tibble::tibble(
				chromosome  = 'Pf3D7_11_v3',
				position    = 1058035,
				start       = 1048035,
				end         = 1068035,
				zoom_start  = 1053000,
				zoom_end    = 1060000
			),
			tree_annotated_positions = c( 1058035, 1057437 )
		),
		`Pfsa3alt` = list(
			tree        = "outputs/pf7/relate/output/trees/popsize/pf7.relate.Pf3D7_11_v3.Ne=100000.mu=6.3162e-9.dpg=365.bp=1057437.newick",
			treesamples = "outputs/pf7/relate/output/trees/popsize/pf7.relate.Pf3D7_11_v3.Ne=100000.mu=6.3162e-9.dpg=365.bp=1057437.samples.newick",
			genotypes   = "outputs/pf7/vcf/07_ancestral/Pf3D7_11_v3.bgen",
			annotation  = "outputs/pf7/vcf/04_merged/Pf3D7_11_v3.merged.annotation.tsv.gz",
			focus       = tibble::tibble(
				chromosome  = 'Pf3D7_11_v3',
				position    = 1057437,
				start       = 1047437,
				end         = 1067437,
				zoom_start  = 1053000,
				zoom_end    = 1060000
			),
			tree_annotated_positions = c( 1057437, 1058035 )
		)
	)
)

data = list()
data$samples = load.samples( args$samples )
data$genes = (
	load.plasmodb.genes( gff = args$genes, gaf = args$gaf )
	%>% filter(
		(ID == 'PF3D7_0215300' | Parent == 'PF3D7_0215300' | Parent == 'PF3D7_0215300.1')
		| (Parent == 'PF3D7_1127000' | ID == 'PF3D7_1127000' | Parent == 'PF3D7_1127000.1')
		| (Parent == 'PF3D7_1126900' | ID == 'PF3D7_1126900' | Parent == 'PF3D7_1126900.1')
		| (ID == 'PF3D7_0220300' | Parent == 'PF3D7_0220300' | Parent == 'PF3D7_0220300.1')
		| (ID == 'PF3D7_0424700' | Parent == 'PF3D7_0424700' | Parent == 'PF3D7_0424700.1')
	)
)

source( "scripts/haplotype_figure_impl.R" )
data$regions = list()
for( name in names( args$regions )) {
	data$regions[[name]] = load.focus.haplotypes.and.tree( args$regions[[name]], data$samples )
}
for( name in names( args$regions )) {
	data$regions[[name]]$plot_spec = create.plot.spec( data$regions[[name]], data$samples, data$genes, subsample.N = 10 )
}

{
	source( "../spatial/code/functions.R" )
	source( "scripts/layout.intervals.R" )
	source( "scripts/plot.genes.R" )
	source( "scripts/haplotype_figure_impl.R" )
	source( "scripts/figure4_impl.R" )
	figure_4(
		specs = list(
			data$regions$Pfsa1$plot_spec,
			data$regions$Pfsa3alt$plot_spec
		),
		colour.column = "Country",
		split = c( 0.425, 0.575 ),
		width = 10,
		height = 10,
		sprintf( "/tmp/figure_4.pdf"),
		colours = list(
			Country = country.colours()[ levels( data$samples$Country ) ],
			region = c(
				west = country.colours()[["Gambia"]],
				central = country.colours()[["Democratic_Republic_of_the_Congo"]],
				east = country.colours()[["Kenya"]]
			)
		)
	)	
}
