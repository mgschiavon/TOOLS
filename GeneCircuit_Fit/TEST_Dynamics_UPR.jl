cd("C:\\Users\\mgsch\\Dropbox (Personal)\\LIIGH\\TOOLS\\GeneCircuit_Fit\\")
using Pkg; Pkg.activate(".");
iARG = (mm = "UPR",	# Label for model file
        ex = "Ex02");	# Label for parameters file

## Load functions & system:
using DelimitedFiles
using DifferentialEquations
using Distributions
mm = include(string("Library\\Md_",iARG.mm,".jl"));	# ODE system
fn = include(string("Library\\FN_Fit.jl"));			# Functions

## INPUTS:
# iARG = (mm : Label for motif file, ex : Label for parameters file);
# System parameters (i.e. p structure -list of kinetic parameters):
include(string("InputFiles\\ARGS_",iARG.mm,"_Par_",iARG.ex,".jl"))
pO = copy(p);

## Load data to compare:
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

##
pOp  = [:gF,:cB,:gUB,:cBI];
vOp  = [0.00595,0.06733,6658.35,0.13356];
for i in 1:length(pOp)
	p[pOp[i]] = vOp[i];
end
p[:cD] = 0.0 * p[:cD0] * p[:mMc] * p[:ERv]; # (DTT [mM])*(cD0*mMc*ERv)
x0d = fn.SS(mm.myODE, p, x0, 1e-6);
plot([0,1000],[x0d[15],x0d[15]],label=string("DTT = 0 mM"),lw=3,legend=:topleft)
for dtt in 2:length(d.DTT)
        p[:cD] = d.DTT[dtt] * p[:cD0] * p[:mMc] * p[:ERv]; # (DTT [mM])*(cD0*mMc*ERv)
        xS = fn.Dyn(mm.myODE, p, x0d, 1000.0, 1e-6);
        plot!(xS.t,[x0d[15],x0d[15]],label=string("DTT = ",d.DTT[dtt]," mM"),lw=3,legend=:topleft)
end
