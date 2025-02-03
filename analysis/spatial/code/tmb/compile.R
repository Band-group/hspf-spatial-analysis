library( argparse )
library( TMB )

echo <- function( message, ... ) {
	cat( sprintf( message, ... ))
}

parse_arguments <- function() {
	parser = ArgumentParser(
		description = 'Compile a TMB model file.'
	)
	parser$add_argument(
		"--model",
		type = "character",
		help = "Path to .cpp file to compile",
		required = TRUE
	)
	return( parser$parse_args() )
}
args = parse_arguments()

echo( "++ Compiling TMB model %s...\n", args$model )

modelfile = gsub( ".cpp", "", args$model )
compile(
	args$model,
	LDFLAGS = "-L/well/band/projects/pfsa-spatial/miniconda/lib/R/lib"	
)
# Compiled ok, but try loading it too...
dyn.load( dynlib( modelfile ))
echo( "++ Compiled." )
