rule fit_hspf_in_areas_with_restricted_sources:
	output:
		rds = "output/pf={pf_data_version}/hspf/fixed-r0={r0}-sigma0={sigma0}-fc={hbs_covariates}/grid-type={type}-size={size}/{locus}/{locus}-model={regression_model}+fc={hspf_covariates}-{min_km_to_survey_pt}km-area={area}-min_N={min_N}-source={source}.rds"
	input:
		grid    = rules.create_grid.output.rds,
		pf      = rules.aggregate_pf.output.tsv,
		hbs     = rules.aggregate_HbS.output.tsv,
		survey  = "input/cleanHbSdata.csv",
		world   = "geodata/naturalearthdata.Rdata"
	params:
		script  = srcdir( "code/BYM-tmb-longform.R" ),
		areas   = lambda w: "" if w.area == 'global' else "--areas '%s'"% "' '".join( config['areas'][w.area] ),
		source  = lambda w: (
			{
				"pf7": ["MalariaGEN Pf7"],
				"moser": ["Moser et al 2021"],
				"verity": ["Verity et al 2021"],
			}[w.source]
		)
	threads: 2
	shell: """
		Rscript --vanilla {params.script} \
		--world {input.world} \
		--grid {input.grid} \
		--model {wildcards.regression_model} \
		--HbS_aggregated {input.hbs} \
		--pf_aggregated {input.pf} \
		--locus {wildcards.locus} \
		{params.areas} \
		--min_km_to_survey_pt {wildcards.min_km_to_survey_pt} \
		--min_N {wildcards.min_N} \
		--sources {params.source} \
		--output {output.rds} \
		--threads {threads}
	"""
