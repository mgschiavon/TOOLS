## ODE FITTING PIPELINE - Gene Regulatory Circuit models
# Fitting instructions - Toy example #1 | Ex01
#	Mariana Gómez-Schiavon
#	August, 2021
#		Julia v.1.5.3
#		Local Environment: "./TOOLS/GeneCircuit_Fit/Project.toml"

using DiffEqBayes
using ParameterizedFunctions, OrdinaryDiffEq, RecursiveArrayTools, Distributions

# Fitting parameters
pOp = Dict([
	:mX => [5.8,5.8*0.2,0,100],       # X synthesis rate (1/min)
	:mY => [2.1,2.1*0.2,0,100],       # Y synthesis rate (1/min)
	:mZ => [0.05,0.05*0.2,0,100],     # Z synthesis rate (1/min)
]);
priors = [];
for i in [eval(Meta.parse(string(":",i))) for i in mm.myODE.sys.ps]
	if(in(i,keys(pOp)))
		priors = vcat(priors,truncated(Normal(pOp[i][1], pOp[i][2]), pOp[i][3], pOp[i][4]));
	else
		priors = vcat(priors,truncated(Normal(p[i], p[i]*1e-4), 0, p[i]*100));
	end
end

# Load data to compare:
# NOTE: d.Xe is always the data matrix used to calculate the MSE.
using CSV
using DataFrames
x = CSV.File("DATA_ToyEx1.csv") |> Tables.matrix;
d = (tD = x[1,2:end],		# Time points
	 Xe = x[2:end,2:end]);	# Data points

# Model:
pV = [p[eval(Meta.parse(string(":",i)))] for i in mm.myODE.sys.ps];
prob1 = ODEProblem(mm.myODE,x0,(0.0,maximum(d.tD)*10),pV);

bayesian_result_hmc = dynamichmc_inference(prob1,Tsit5(),d.tD,d.Xe,priors)

# Only Z:
x = CSV.File("DATA_ToyEx1.csv") |> Tables.matrix;
d = (tD = x[1,2:end],		# Time points
	 Xe = x[end,2:end]);	# Data points

bayesian_result_hmc2 = dynamichmc_inference(prob1,Tsit5(),d.tD,reshape(d.Xe,1,length(d.Xe)),priors,save_idxs=[3])
