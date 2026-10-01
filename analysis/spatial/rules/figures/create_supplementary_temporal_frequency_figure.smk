rule create_supplementary_temporal_frequency_figure:
	output:
		pdf = "output/pf={pf_data_version}/figures/temporal/{loci}-temporal-area={area}.pdf"
	input:
		tsv = "output/pf={pf_data_version}/pf/aggregated/grid-type=hexagon-size=1-area={area}-by=year-source.tsv"
	params:
		script = srcdir( "code/figures/temporal_figure.R" ),
		loci = lambda w: w.loci.split( "+" ),
		#countries = lambda w: (
		#	config['areas'].get( w.area, w.area.replace( "+", " " ) )
		#)
		# We have hard-coded this for the figure.
		countries = "--countries %s" % ' '.join([
			"Gambia",
			"Senegal",
			"Mali",
			"Ghana",
			"Nigeria",
			"Democratic_Republic_of_the_Congo",
			"Uganda",
			"Malawi",
			"Mozambique",
			"Kenya",
			"Tanzania"
		])
	shell: """
	Rscript --vanilla {params.script} \
	--pf_aggregated {input.tsv} \
	--loci {params.loci} \
	{params.countries} \
	--output {output.pdf}
"""
