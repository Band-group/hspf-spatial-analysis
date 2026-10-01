#' Load genotypes and read counts for specified variants from a VCF file
#'
#' @description
#' Load genotypes and read counts for specified variants from a VCF file,
#' Uses bcftools query under the hood; this must be installed.
#'
#' @param filename Path to an input file.
#' @param variants A dataframe of variants with chromosome and position columns
#'
#' @return a data frame of samples, genotypes, and ref and alt read counts
#'
#' @export
load.genotypes.from.vcf <- function( filename, variants ) {
	library( dplyr )
	t1 = tempfile()
	data = readr::read_tsv( 
		pipe(
			sprintf(
				"bcftools query -f '[%%SAMPLE\t%%CHROM\t%%POS\t%%REF\t%%ALT\t%%GT\t%%AD\n]' %s",
				filename
			)
		),
		col_names = c( "ID", "chromosome", "position", "ref", "alt", "GT", "AD"),
		col_types = c( "ccicccc" )
	)
	stopifnot( nrow( data ) > 0 )
	data$GT[ data$GT == './.' ] = NA

	# Fix chromosomes written as chr%d.
	data$chromosome_number = as.integer( stringr::str_extract( data$chromosome, '(Pf3D7_|chr)([0-9]+)(_v3|)', group = 2 ) )
	data$chromosome = sprintf( "Pf3D7_%02d_v3", data$chromosome_number )

	data = (
		data
		%>% inner_join(
			variants,
			by = c( 'chromosome', 'position' )
		)
	)
	data = (
		data
		%>% mutate(
			read_count_ref = as.integer( stringr::str_extract( AD, "([0-9]+),([0-9]+)(,[0-9]+|)", group = 1 )),
			read_count_alt = as.integer( stringr::str_extract( AD, "([0-9]+),([0-9]+)(,[0-9]+|)", group = 2 ))
		)
	)
	return(data)
}
