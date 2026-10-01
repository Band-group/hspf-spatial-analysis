# We compile the C++ code to a directory named by the
# current system's machine (arm64 or x86_64).
# This avoids trying to use the wrong .so file if we copy the output directory.
import platform

rule compile_TMB_code:
	output:
		cpp = "output/hspf/tmb/{platform}/{regression_model}.cpp".format( platform = platform.machine(), regression_model = '{regression_model}' ),
		so  = "output/hspf/tmb/{platform}/{regression_model}.so".format( platform = platform.machine(), regression_model = '{regression_model}' ),
	input:
		cpp = srcdir( "code/tmb/{regression_model}.cpp" )
	params:
		libpath = "/well/band/projects/pfsa-spatial/miniconda/lib/R/lib/",
		script = srcdir( "code/tmb/compile.R" )
	shell: """
		cp {input.cpp} {output.cpp}
		Rscript --vanilla {params.script} --model {output.cpp}
	"""
