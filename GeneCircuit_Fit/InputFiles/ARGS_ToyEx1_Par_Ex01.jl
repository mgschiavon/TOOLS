## ODE FITTING PIPELINE - Gene Regulatory Circuit models
# Parameters & other instructions - Toy example #1 | Ex01
#	Mariana Gómez-Schiavon
#	August, 2021
#		Julia v.1.5.3
#		Local Environment: "./TOOLS/GeneCircuit_Fit/Project.toml"

# Kinetic parameters
p = Dict([
	:mX => 1.0,       # X synthesis rate (1/min)
	:gX => 0.01,      # X degradation rate (1/min)
	:mY => 0.2,       # Y synthesis rate (1/min)
	:gY => 0.01,      # Y degradation rate (1/min)
	:mZ => 0.15,      # Z synthesis rate (1/min)
	:gZ => 0.01,      # Z degradation rate (1/min)
]);

# Initial conditions:
x0 = zeros(length(mm.myODE.syms));
x0[1] = 0;		  	  # X
x0[2] = 0;           # Y
x0[3] = 0;           # Z
