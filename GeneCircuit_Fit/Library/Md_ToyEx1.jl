## ODE FITTING PIPELINE - Gene Regulatory Circuit models
# ODE model - Toy example #1
#	Mariana Gómez-Schiavon
#	August, 2021
#		Julia v.1.5.3
#		Local Environment: "./TOOLS/GeneCircuit_Fit/Project.toml"

module mm
	# Required libraries
	using DifferentialEquations
	using ParameterizedFunctions

	# ODE system
	myODE = @ode_def begin
		dX  = mX - (gX * X)
		dY  = (mY * X) - (gY * Y)
		dZ  = (mZ * Y) - (gZ * Z)
	end mX gX mY gY mZ gZ;
end
