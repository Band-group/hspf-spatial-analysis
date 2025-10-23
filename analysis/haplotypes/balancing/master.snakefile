
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

countries = [
	"Gambia",
	"Senegal",
#	"Guinea",
#	"Mauritania",
#	"Cote_dIvoire",
	"Mali",
#	"Burkina_Faso",
	"Ghana",
	"Benin",
#	"Nigeria",
#	"Gabon",
	"Cameroon",
	"Democratic_Republic_of_the_Congo",
#	"Sudan",
#	"Uganda",
	"Malawi",
	"Tanzania",
#	"Mozambique",
	"Kenya",
#	"Ethiopia",
#	"Madagascar
]

country_sets = {
	'east': [ 'Kenya', 'Malawi', 'Tanzania', 'Uganda', 'Mozambique' ],
	'west': [ 'Gambia', 'Senegal', 'Guinea', 'Mauritania', 'Cote_dIvoire', 'Burkina_Faso', 'Ghana', 'Benin', 'Mali' ],
	'central_west': [ 'Nigeria', 'Cameroon', 'Gabon', 'Democratic_Republic_of_the_Congo' ],
	'north_east': [ 'Sudan', 'Ethiopia' ]
}

year_sets = {
	'1980-1999': { 'start': 1960, 'end': 1999 },
	'2000-2009': { 'start': 2000, 'end': 2009 },
	'2010-2025': { 'start': 2010, 'end': 2026 },
	'all': { 'start': 1960, 'end': 2026 }
}

wildcard_constraints:
	chromosome = "|".join( chromosomes ),
	Ne = '[0-9]+',
	beagle_version = "v5[.]4|v5[.]1|v4[.]1"

include: "rules/pf7.snakefile"
include: "rules/selscan.snakefile"
include: "rules/betascan.snakefile"

localrules: combine_selscan

rule all:
	input:
		vcf = expand( "outputs/pf7/vcf/06_phased/{chromosome}.phased.v5.4.vcf.gz", chromosome = chromosomes ),
		counts = expand( "outputs/pf7/vcf/06_phased/{chromosome}.phased.{beagle_version}.counts.txt", chromosome = chromosomes, beagle_version = [ "v5.4" ]),
		ancestral = expand( "outputs/pf7/vcf/07_ancestral/{chromosome}.{extension}", chromosome = chromosomes, extension = [ 'bgen', 'vcf.gz', 'shapeit.gz' ] ),
		#polarised = expand( "outputs/pf7/relate/input/{chromosome}.shapeit.gz", chromosome = chromosomes ),
		#samples = "outputs/pf7/relate/input/relate_input.sample",
		betascan = expand(
			"outputs/pf7/betascan/output/pf7.betascan.window={window}.p={p}.tsv.gz",
			window = [ "5000", "10000" ],
			p = [ "20", "50" ]
		),
		selscan = expand(
			"outputs/pf7/selscan/output/pf7.selscan.{mode}.bins={bins}.tsv.gz",
			mode = [ 'ihs', 'ihh12' ],
			bins = [ '1%', '2.5%', '5%' ]
		)
