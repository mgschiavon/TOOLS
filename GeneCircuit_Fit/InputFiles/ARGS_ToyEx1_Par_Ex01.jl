## ODE FITTING PIPELINE - Gene Regulatory Circuit models
# Parameters & fitting instructions - Toy example #1 | Ex01
#	Mariana Gómez-Schiavon
#	August, 2021
#		Julia v.1.5.3
#		Local Environment: "./TOOLS/GeneCircuit_Fit/Project.toml"

# Kinetic parameters
p = Dict([
	:mX => 5.8,       # X synthesis rate (1/min)
	:gX => 0.01,      # X degradation rate (1/min)
	:mY => 2.1,       # Y synthesis rate (1/min)
	:gY => 0.01,      # Y degradation rate (1/min)
	:mZ => 0.05,      # Z synthesis rate (1/min)
	:gZ => 0.01,      # Z degradation rate (1/min)
]);

# Initial conditions:
x0 = zeros(3);
x0[1] = 0;		  	 # X
x0[2] = 0;           # Y
x0[3] = 0;           # Z

# Fitting instructions ([mean, standard deviation, minimum value, maximum value]):
pOp = Dict([
	:mX => [5.8,5.8*0.2,0,100],       # X synthesis rate (1/min)
	:mY => [2.1,2.1*0.2,0,100],       # Y synthesis rate (1/min)
	:mZ => [0.05,0.05*0.2,0,100],     # Z synthesis rate (1/min)
]);

# Reference data:
x  = CSV.File("DATA_ToyEx1.csv") |> Tables.matrix;
tD = x[1,2:end];		# Time points
vD = [1,2,3];			# Observable variables
xD = x[2:end,2:end];	# Data points (array size: vD x tD)
