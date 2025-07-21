chromosomes = ["Pf3D7_%02d_v3" % i for i in range( 1, 15 )]

data = {
	"pf7": "/well/band/projects/pf7"
}

regions = {
	"Pfsa1": {
		"foci": [ "631190" ],
		"chromosome": "Pf3D7_02_v3",
		"range": "Pf3D7_02_v3:581190-681190",
		"exclusions": [
			'Pf3D7_02_v3:630737' # This variant appears to be a gene conversion or similar
		]
	},
	"Pfsa2": {
		"foci": [ "814288" ],
		"chromosome": "Pf3D7_02_v3",
		"range": "Pf3D7_02_v3:764288-864288",
		"exclusions": []
	},
	"Pfsa3": {
		"foci": [ "1058035", "1057437" ],
		"chromosome": "Pf3D7_11_v3",
		"range": "Pf3D7_11_v3:1008035-1108035",
		"exclusions": [
			# R2
			"Pf3D7_11_v3:1054003",
			"Pf3D7_11_v3:1054020",
			"Pf3D7_11_v3:1054021",
			"Pf3D7_11_v3:1054078",
			"Pf3D7_11_v3:1054099",
			"Pf3D7_11_v3:1054106",
			"Pf3D7_11_v3:1054229",
			"Pf3D7_11_v3:1054232",
			"Pf3D7_11_v3:1054282",
			"Pf3D7_11_v3:1054287",
			"Pf3D7_11_v3:1054303",
			"Pf3D7_11_v3:1054310",
			"Pf3D7_11_v3:1054343",
			"Pf3D7_11_v3:1054444",
			"Pf3D7_11_v3:1054447",
			"Pf3D7_11_v3:1054587",
			"Pf3D7_11_v3:1054608",
			"Pf3D7_11_v3:1054633",
			"Pf3D7_11_v3:1054634",
			"Pf3D7_11_v3:1054655",
			"Pf3D7_11_v3:1054659",
			"Pf3D7_11_v3:1054683",
			"Pf3D7_11_v3:1054711",
			"Pf3D7_11_v3:1054876",
			"Pf3D7_11_v3:1054880",
			"Pf3D7_11_v3:1054882",
			"Pf3D7_11_v3:1054944",
			"Pf3D7_11_v3:1055033",
			# R3
			"Pf3D7_11_v3:1055085",
			"Pf3D7_11_v3:1055103",
			"Pf3D7_11_v3:1055104",
			"Pf3D7_11_v3:1055108",
			"Pf3D7_11_v3:1055152",
			"Pf3D7_11_v3:1055223",
			"Pf3D7_11_v3:1055225",
			"Pf3D7_11_v3:1055275",
			"Pf3D7_11_v3:1055379",
			"Pf3D7_11_v3:1055408",
			# R4
			"Pf3D7_11_v3:1055455",
			"Pf3D7_11_v3:1055693",
			# R7
			"Pf3D7_11_v3:1058826",
			"Pf3D7_11_v3:1058828",
			"Pf3D7_11_v3:1058874",
			"Pf3D7_11_v3:1058895",
			"Pf3D7_11_v3:1058971",
			"Pf3D7_11_v3:1058996",
			"Pf3D7_11_v3:1059013",
			# R8
			"Pf3D7_11_v3:1059146",
			"Pf3D7_11_v3:1059172",
			"Pf3D7_11_v3:1059173",
			"Pf3D7_11_v3:1059187",
			"Pf3D7_11_v3:1059193",
			"Pf3D7_11_v3:1059200",
			"Pf3D7_11_v3:1059205",
			"Pf3D7_11_v3:1059315",
			"Pf3D7_11_v3:1059357",
			"Pf3D7_11_v3:1059399",
			"Pf3D7_11_v3:1059446",
			"Pf3D7_11_v3:1059448",
			"Pf3D7_11_v3:1059464",
			"Pf3D7_11_v3:1059465",
			"Pf3D7_11_v3:1059477",
			"Pf3D7_11_v3:1059491",
			"Pf3D7_11_v3:1059494",
			"Pf3D7_11_v3:1059510",
			"Pf3D7_11_v3:1059553",
			"Pf3D7_11_v3:1059560",
			"Pf3D7_11_v3:1059577",
			"Pf3D7_11_v3:1059608",
			"Pf3D7_11_v3:1059621",
			"Pf3D7_11_v3:1059635",
			# R9
			"Pf3D7_11_v3:1059737"
			# R10
		]
	}
}
wildcard_constraints:
	chromosome = "|".join( chromosomes ),
	position = '[0-9]+',
	Ne = '[0-9]+'

include: "relate.snakefile"

mutation_rates = [
	# Otto et al says average mutation rate is 9.57x10-11 per mitosis
	# (which comes from the average 3.83x10^-10 per erythrocyte cycle from Claessens et al clone tree paper, and 4 mitoses per erythrocyte cycle.)
	# The table then says between 66 and 336 mitosis (i.e. about 16-84 erythocyte cycles) per generation
	# so the mutation rate per site per generation this is:
	'6.3162e-9', '3.21552e-8'
	# We previously used 4.35e-9 which is evidently on the low side.
]

rule all:
	input:
		relate = expand(
			"outputs/pf7/relate/output/initial/pf7.relate.{chromosome}.Ne=100000.mu={mu}.mut",
			chromosome = chromosomes,
			# Otto et al says average mutation rate is 9.57x10-11 per mitosis
			# (which comes from the average 3.83x10^-10 per erythrocyte cycle from Claessens et al clone tree paper, and 4 mitoses per erythrocyte cycle.)
			# The table then says between 66 and 336 mitosis (i.e. about 16-84 erythocyte cycles) per generation
			# so the mutation rate per site per generation this is:
			mu = mutation_rates
			# We previously used 4.35e-9 which is evidently on the low side.
		),
		#regions = expand( "outputs/pf7/relate/input/{region}.shapeit.gz", region = regions.keys() ),
		popsize = expand(
			"outputs/pf7/relate/output/popsize/pf7.relate.{chromosome_or_region}.Ne={Ne}.mu={mu}.popsize.pdf",
			chromosome_or_region = chromosomes,
			mu = mutation_rates,
			Ne = [ "100000" ]
		)
