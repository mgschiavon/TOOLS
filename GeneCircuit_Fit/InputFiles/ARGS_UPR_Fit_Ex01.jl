## ODE FITTING PIPELINE - Gene Regulatory Circuit models
# Fitting instructions - Unfolded protein response | Ex01
#	Mariana Gómez-Schiavon
#	August, 2021
#		Julia v.1.5.3
#		Local Environment: "./TOOLS/GeneCircuit_Fit/Project.toml"

# Fitting parameters
mrw  = (pOp  = [:gF,:cB,:gUB,:cBI],# Parameters to optimize
		pMin = [-5,-4,0,-4],# Minimum parameter value to explore (log10)
		pMax = [-1,0,4,0],# Maximum parameter value to explore (log10)
		runs = 4,# Number of optimization runs
		iter = 100000,# Number of iterations per optimization run
		cov  = [0.1,0.1,0.1],# Covariance to calculate parameter random walk
		M    = 2,	# "Mutation step size" for multiplicative random walk
		rnP0 = 0,	# Flag for random initial values of parameters to optimize
		temp = 0,	# Flag for simulated annealing (if 0, MRW)
		prtW = 1);	# Flag for printing each walk step

# Load data to compare:
# NOTE: d.Xe is always the data matrix used to calculate the MSE.
using CSV
using DataFrames
### Using estimated data (from Pincus et al., 2010, Fig. 3) on           ###
###    '\MODEL - Quantifying feedback\RUN_Natural_UPR\Data_Fig3.xlsx',   ###
###    with rows 1,9,17 for (same) time points, and rows 2-8 for Fig. 3A,###
###    rows 10-16 for Fig. 3B, and rows 18-24 for Fig. 3C data points.   ###
x = CSV.File("DATA_UPR.csv") |> Tables.matrix;
d = (tD  = x[1,2:end],						# Time points (min)
	 DTT = parse.(Float64, x[2:8,1]),		# DTT concetration (mM)
	 Xe  = x[vcat(2:8,10:16,18:24),2:end]);	# Data points

# RULES:
function mySIM(fn,mm,p,d)
	Y = zeros(size(d.Xe));
	pWT = copy(p);
	# Calculate dynamics -- WT (Fig. 3A, Pincus et al.):
	for dtt in 1:length(d.DTT)
		p[:cD] = d.DTT[dtt] * p[:cD0] * p[:mMc] * p[:ERv]; # (DTT [mM])*(cD0*mMc*ERv)
		xS = fn.Dyn(mm.myODE, p, x0, 300.0, 1e-6);
		for i in 1:length(d.tD)
			Y[dtt,i] = xS(d.tD[i])[15];
		end
	end
	# Calculate dynamics -- hac1 mutant (Fig. 3B, Pincus et al.):
	p = copy(pWT);
	p[:nB] = 0;
	p[:nE] = 0;
	for dtt in 1:length(d.DTT)
		p[:cD] = d.DTT[dtt] * p[:cD0] * p[:mMc] * p[:ERv]; # (DTT [mM])*(cD0*mMc*ERv)
		xS = fn.Dyn(mm.myODE, p, x0, 300.0, 1e-6);
		for i in 1:length(d.tD)
			Y[dtt+length(d.DTT),i] = xS(d.tD[i])[15];
		end
	end
	# Calculate dynamics -- Ire1-bipless (Fig. 3C, Pincus et al.):
	p = copy(pWT);
	p[:cBI] = 0;
	p[:gIB] = 0;
	for dtt in 1:length(d.DTT)
		p[:cD] = d.DTT[dtt] * p[:cD0] * p[:mMc] * p[:ERv]; # (DTT [mM])*(cD0*mMc*ERv)
		xS = fn.Dyn(mm.myODE, p, x0, 300.0, 1e-6);
		for i in 1:length(d.tD)
			Y[dtt+(2*length(d.DTT)),i] = xS(d.tD[i])[15];
		end
	end
	p = copy(pWT);
	return Y;
end
