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
# Fitting rules (i.e. mrw structure -list of metropolis random walk parameters;
#					  d structure -experimental data/conditions;
#					  mySIM -specific function to calculate simulated data):
include(string("InputFiles\\ARGS_",iARG.mm,"_Fit_",iARG.ex,".jl"))
pO = copy(p);

## Test dynamics:
pOp  = [:gF,:cB,:gUB,:cBI];
vOp  = [0.00595,0.06733,6658.35,0.13356];
for i in 1:length(pOp)
	p[pOp[i]] = vOp[i];
end
Y = mySIM(fn,mm,p,d)
open(string("TEST_DYN.txt"), "w") do io
	writedlm(io, [vcat("DTT",[string(t) for t in d.tD])],'\t');
	for i in 1:size(Y)[1]
		writedlm(io, [vcat(vcat(d.DTT,d.DTT,d.DTT)[i],Y[i,:])],'\t');
	end
end

## Test stability:
d.tD .*= 10;
Y = mySIM(fn,mm,p,d)
open(string("TEST_DYN_LONG.txt"), "w") do io
	writedlm(io, [vcat("DTT",[string(t) for t in d.tD])],'\t');
	for i in 1:size(Y)[1]
		writedlm(io, [vcat(vcat(d.DTT,d.DTT,d.DTT)[i],Y[i,:])],'\t');
	end
end
