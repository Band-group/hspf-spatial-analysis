rule hspf_summaries_to_excel:
	output:
		xlsx = rules.combine_hspf_summaries.output.tsv.replace( '.tsv', '.xlsx' )
	input:
		tsv = rules.combine_hspf_summaries.output.tsv
	params:
		script = srcdir( "code/summary_hspf2excel.R" )
	shell: """
		Rscript --vanilla {params.script} \
		--input {input.tsv} \
		--output {output.xlsx}
	"""
