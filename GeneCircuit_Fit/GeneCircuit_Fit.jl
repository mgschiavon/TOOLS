## ODE FITTING PIPELINE - Gene Regulatory Circuit models
#  Main pipeline
#	Mariana Gómez-Schiavon
#	August, 2021
#		Julia v.1.5.3
#		Local Environment: "./TOOLS/GeneCircuit_Fit/Project.toml"

## Load functions & system:
using RecursiveArrayTools
using Distributions
using DelimitedFiles
using CSV
using DataFrames
using DifferentialEquations
using ParameterizedFunctions
using DiffEqBayes
using OrdinaryDiffEq
using DynamicHMC

## INPUTS:
# iARG = (mm : Label for motif file, ex : Label for parameters file);
# Parameters & fitting instructions:
include(string("InputFiles\\ARGS_",iARG.mm,"_Par_",iARG.ex,".jl"))

## Define ODE system
par  = Dict([
	:Fix => "";		# Fixed parameters
	:Fit => "end ";	# Parameters to be fitted
	:priors => [];		# Prior distributions
	]);
for i in keys(p)
	if(in(i,keys(pOp)))
		par[:Fit] *= String(i)*" ";
		par[:priors]  = vcat(par[:priors],truncated(Normal(pOp[i][1], pOp[i][2]), pOp[i][3], pOp[i][4]));
	else
		par[:Fix] *= String(i)*"="*string(p[i])*"; ";
	end
end

eqs = readlines(string("Library\\Md_",iARG.mm,".jl"));
open("InputFiles\\myODE_Syst.jl","w") do io
	writedlm(io, [vcat(par[:Fix])],'\n')
	writedlm(io, [vcat("myODE = @ode_def begin")],'\n');
	writedlm(io, [vcat(eqs)],'\n')
	writedlm(io, [vcat(par[:Fit])],'\n')
end
include("InputFiles\\myODE_Syst.jl")
rm("InputFiles\\myODE_Syst.jl")

## Dynamic Hamilton Monte Carlo fitting:
pV = [p[eval(Meta.parse(string(":",i)))] for i in myODE.sys.ps];
hm = dynamichmc_inference(ODEProblem(myODE,x0,(0.0,maximum(tD)*10.0),pV),Tsit5(),tD,xD,par[:priors],save_idxs=vD);
# Print output:
open(string("OUT_Fit_",iARG.mm,"_",iARG.ex,".cvs"), "w") do io
	writedlm(io, [vcat([i for i in keys(pOp)],["s_"*string(myODE.syms[i]) for i in vD])],'\t');
	for i in 1:length(hm.posterior)
		writedlm(io, [vcat(hm.posterior[i].parameters,hm.posterior[i].σ)],'\t')
	end
end
