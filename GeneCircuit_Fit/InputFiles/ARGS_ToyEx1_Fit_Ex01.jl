## ODE FITTING PIPELINE - Gene Regulatory Circuit models
# Fitting instructions - Toy example #1 | Ex01
#	Mariana Gómez-Schiavon
#	August, 2021
#		Julia v.1.5.3
#		Local Environment: "./TOOLS/GeneCircuit_Fit/Project.toml"

# Fitting parameters
mrw  = (pOp  = [:mX,:mY,:mZ],# Parameters to optimize
		pMin = [-2,-2,-2],# Minimum parameter value to explore (log10)
		pMax = [2,2,2],# Maximum parameter value to explore (log10)
		runs = 4,# Number of optimization runs
		iter = 100000,# Number of iterations per optimization run
		cov  = [0.01,0.01,0.01],# Covariance to calculate parameter random walk
		M    = 2,	# "Mutation step size" for multiplicative random walk
		rnP0 = 0,	# Flag for random initial values of parameters to optimize
		temp = 0,	# Flag for simulated annealing (if 0, MRW)
		prtW = 1);	# Flag for printing each walk step

# Load data to compare:
# NOTE: d.Xe is always the data matrix used to calculate the MSE.
using CSV
using DataFrames
x = CSV.File("DATA_ToyEx1.csv") |> Tables.matrix;
d = (tD = x[1,2:end],		# Time points
	 Xe = x[2:end,2:end]);	# Data points

# RULES:
function mySIM(fn,mm,p,d)
	Y = zeros(size(d.Xe));
	# Calculate dynamics:
	xS = fn.Dyn(mm.myODE, p, x0, 10.0, 1e-6);
	for i in 1:length(d.tD)
		Y[:,i] = xS(d.tD[i]);
	end
	return Y;
end
