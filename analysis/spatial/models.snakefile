ranges = [
	'5', '10', '15'
]
sigmas = [
	'0.6', '0.8', '1.0'
]
covariates = [ 'none', 'continent' ]

rule all:
	input:
		filenames = expand(
			"output/HbSsensitivity/fits/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}/filenames.tsv",
			r0 = ranges,
			sigma0 = sigmas,
			covariates = covariates
		)

rule fit_hbs_map:
	output:
		filenames	= "output/HbSsensitivity/fits/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}/filenames.tsv",
		fit 		= "output/HbSsensitivity/fits/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}_modelfit.rds",
		predictions	= "output/HbSsensitivity/fits/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}_predictions.rds",
		prior		= "output/HbSsensitivity/fits/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}_prior.tsv",
		samples		= "output/HbSsensitivity/fits/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}_samples.rds",
		xyt			= "output/HbSsensitivity/fits/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}_xyt.rds"
	input:
		hbs = "input/cleanHbSdata.csv",
		piel = 'geodata/2013_Sickle_Haemoglobin_HbS_Allele_Freq_Global_5k_Decompressed.tif',
		geodata = directory('geodata')
	params:
		script = "code/HbS_model_fit2.R",
		outdir = "output/HbSsensitivity/fits/fixed-r0={r0}-sigma0={sigma0}-fc={covariates}",
		covariates = lambda wildcards: ( '--fixed_covariates %s' % wildcards.covariates if wildcards.covariates != 'none' else '' )
	shell: """
	Rscript --vanilla {params.script} \
	--geodata geodata \
	--HbS input/cleanHbSdata.csv \
	--piel geodata/2013_Sickle_Haemoglobin_HbS_Allele_Freq_Global_5k_Decompressed.tif \
	--r0 {wildcards.r0} \
	--sigma0 {wildcards.sigma0} \
	{params.covariates} \
	--outdir {params.outdir}
"""
