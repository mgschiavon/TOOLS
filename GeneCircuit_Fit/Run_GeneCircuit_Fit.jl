## Running in julia terminal
	cd("\\TOOLS\\GeneCircuit_Fit\\")	# Go to files location
	using Pkg; Pkg.activate(".");		# Activate local environment (requiere '.toml' files)
	iARG = (mm = "ToyEx1",				# Label for model file (e.g. 'A' for 'Library\Md_A.jl' file)
			ex = "Ex01");				# Label for parameters file (e.g. 'Ex01' for 'InputFiles\ARGS_A_Par_Ex01.jl' file)
	include("GeneCircuit_Fit.jl")		# Run fiting code
