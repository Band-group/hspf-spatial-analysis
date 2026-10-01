rule fit_hbs_map:
	output:
		filenames	= "output/HbS/{hbs_model_type}-r0={r0}-sigma0={sigma0}-fc={hbs_covariates}/fit/catalogue.tsv",
		prior		= "output/HbS/{hbs_model_type}-r0={r0}-sigma0={sigma0}-fc={hbs_covariates}/fit/{hbs_model_type}-r0={r0}-sigma0={sigma0}-fc={hbs_covariates}_prior.tsv",
		xyt			= "output/HbS/{hbs_model_type}-r0={r0}-sigma0={sigma0}-fc={hbs_covariates}/fit/{hbs_model_type}-r0={r0}-sigma0={sigma0}-fc={hbs_covariates}_xyt.rds",
		fit 		= "output/HbS/{hbs_model_type}-r0={r0}-sigma0={sigma0}-fc={hbs_covariates}/fit/{hbs_model_type}-r0={r0}-sigma0={sigma0}-fc={hbs_covariates}_modelfit.rds",
		predictions	= "output/HbS/{hbs_model_type}-r0={r0}-sigma0={sigma0}-fc={hbs_covariates}/fit/{hbs_model_type}-r0={r0}-sigma0={sigma0}-fc={hbs_covariates}_predictions.rds",
		samples		= "output/HbS/{hbs_model_type}-r0={r0}-sigma0={sigma0}-fc={hbs_covariates}/fit/{hbs_model_type}-r0={r0}-sigma0={sigma0}-fc={hbs_covariates}_samples.rds",
#		HbSmesh		= "output/HbS/{hbs_model_type}-r0={r0}-sigma0={sigma0}-fc={hbs_covariates}/fit/HbSmesh.pdf"
	input:
		hbs = "input/cleanHbSdata.csv",
		piel = 'geodata/2013_Sickle_Haemoglobin_HbS_Allele_Freq_Global_5k_Decompressed.tif',
		geodata = directory('geodata')
	params:
		script         = "code/HbS_model_fit2.R",
		outdir         = "output/HbS/{hbs_model_type}-r0={r0}-sigma0={sigma0}-fc={hbs_covariates}/fit",
		hbs_covariates = lambda wildcards: ( '--fixed_covariates %s' % wildcards.hbs_covariates if wildcards.hbs_covariates != 'none' else '' ),
		model_options  = lambda w: ['--model_options +iid'] if '+iid' in w.hbs_model_type else ''
	wildcard_constraints:
		hbs_model_type = 'fixed|fixed[+]iid',
		r0             = '[^()]+',
		sigma0         = '[^()]+'
	shell: """
	Rscript --vanilla {params.script} \
	--geodata geodata \
	--HbS input/cleanHbSdata.csv \
	--piel geodata/2013_Sickle_Haemoglobin_HbS_Allele_Freq_Global_5k_Decompressed.tif \
	--r0 {wildcards.r0} \
	--sigma0 {wildcards.sigma0} \
	{params.model_options} \
	{params.hbs_covariates} \
	--outdir {params.outdir}
"""

# FOr PC prior specify r0 and sigma0 as:
# r0=PC(r,P) where r is range the parameter and P is the PC prior probability
# sigma0=PC(sigma,P) where sigma is the variance parameter and P is the PC prior probability
# Note make sure to use quotes around filenames due to the brackets.
rule fit_hbs_map_with_PC_prior:
	output:
		filenames	   = "output/HbS/{hbs_model_type}-r0={r0}-sigma0={sigma0}-fc={hbs_covariates}/fit/catalogue.tsv",
		prior		   = "output/HbS/{hbs_model_type}-r0={r0}-sigma0={sigma0}-fc={hbs_covariates}/fit/{hbs_model_type}-r0={r0}-sigma0={sigma0}-fc={hbs_covariates}_prior.tsv",
		xyt			   = "output/HbS/{hbs_model_type}-r0={r0}-sigma0={sigma0}-fc={hbs_covariates}/fit/{hbs_model_type}-r0={r0}-sigma0={sigma0}-fc={hbs_covariates}_xyt.rds",
		fit 		   = "output/HbS/{hbs_model_type}-r0={r0}-sigma0={sigma0}-fc={hbs_covariates}/fit/{hbs_model_type}-r0={r0}-sigma0={sigma0}-fc={hbs_covariates}_modelfit.rds",
		predictions	   = "output/HbS/{hbs_model_type}-r0={r0}-sigma0={sigma0}-fc={hbs_covariates}/fit/{hbs_model_type}-r0={r0}-sigma0={sigma0}-fc={hbs_covariates}_predictions.rds",
		samples		   = "output/HbS/{hbs_model_type}-r0={r0}-sigma0={sigma0}-fc={hbs_covariates}/fit/{hbs_model_type}-r0={r0}-sigma0={sigma0}-fc={hbs_covariates}_samples.rds",
#		HbSmesh		   = "output/HbS/fixed-r0={r0}-sigma0={sigma0}-fc={hbs_covariates}/fit/HbSmesh.pdf"
	input:
		hbs            = "input/cleanHbSdata.csv",
		piel           = 'geodata/2013_Sickle_Haemoglobin_HbS_Allele_Freq_Global_5k_Decompressed.tif',
		geodata        = directory('geodata')
	params:
		script         = "code/HbS_model_fit2.R",
		outdir         = "output/HbS/{hbs_model_type}-r0={r0}-sigma0={sigma0}-fc={hbs_covariates}/fit",
		hbs_covariates = lambda wildcards: ( '--fixed_covariates %s' % wildcards.hbs_covariates if wildcards.hbs_covariates != 'none' else '' ),
		r0_spec        = lambda w: w.r0.replace( 'PC(', '' ).replace( ')', '' ).split( ',' ),
		sigma0_spec    = lambda w: w.sigma0.replace( 'PC(', '' ).replace( ')', '' ).split( ',' ),
		model_options  = lambda w: ['--model_options +iid'] if '+iid' in w.hbs_model_type else ''
	wildcard_constraints:
		hbs_model_type = 'variable|variable[+]iid',
		r0     = 'PC[(].*,.*[)]',
		sigma0 = 'PC[(].*,.*[)]'
	shell: """
	echo {params.r0_spec}
	echo {params.sigma0_spec}
	Rscript --vanilla {params.script} \
	--geodata geodata \
	--HbS input/cleanHbSdata.csv \
	--piel geodata/2013_Sickle_Haemoglobin_HbS_Allele_Freq_Global_5k_Decompressed.tif \
	--r0 {params.r0_spec[0]} \
	--Prange {params.r0_spec[1]} \
	--sigma0 {params.sigma0_spec[0]} \
	--Psigma {params.sigma0_spec[1]} \
	{params.model_options} \
	{params.hbs_covariates} \
	--outdir '{params.outdir}'
"""
