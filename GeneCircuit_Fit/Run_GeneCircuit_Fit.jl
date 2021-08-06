## Running in julia terminal
	cd("C:\\Users\\mgsch\\Dropbox (Personal)\\LIIGH\\TOOLS\\GeneCircuit_Fit\\")
	using Pkg; Pkg.activate(".");
	iARG = (mm = "UPR",	# Label for model file
			ex = "Ex01");	# Label for parameters file
	include("GeneCircuit_Fit.jl")
