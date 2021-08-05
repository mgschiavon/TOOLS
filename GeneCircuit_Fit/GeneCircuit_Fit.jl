## ODE FITTING PIPELINE - Gene Regulatory Circuit models
#  Main pipeline
#	Mariana Gómez-Schiavon
#	August, 2021
#		Julia v.1.5.3
#		Local Environment: "./TOOLS/GeneCircuit_Fit/Project.toml"

# Load functions & system:
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

# Run analysis
open(string("OUT_Fit_",iARG.mm,"_",iARG.ex,".txt"), "w") do io
	writedlm(io, [vcat("Run","Iteration","MSE",[string(param) for param in mrw.pOp])],'\t');
	for ruN in 1:mrw.runs
		println("RUN #",ruN)
		## Reset parameters to input values:
		p = copy(pO);
		## Random initial values of parameters to optimize:
		if(mrw.rnP0==1)
			for i in 1:length(mrw.pOp)
				p[mrw.pOp[i]] = 10 .^ (rand(Uniform(mrw.pMin[i], mrw.pMax[i])));
			end
		end
		## Temperature function for simulated annealing:
		if(mrw.temp==1)
			mrwT = collect(mrw.iter:-1:1) ./ mrw.iter;
		else
			mrwT = ones(mrw.iter); # NOTE: For MRW, make T=1.
		end
		## Initialize system:
		# First calculation of simulated data + MSE:
		mse0 = fn.MSE(mySIM(fn,mm,p,d),d.Xe);
		r0 = zeros(length(mrw.pOp));
		writedlm(io, [vcat(ruN,0,mse0,[p[i] for i in mrw.pOp])],'\t')
		# Optimization iterations:
		println("I: MSE = ",mse0)
		for i in 1:mrw.iter
			# Random values to update parameters:
			rI = rand(MvNormal(zeros(length(mrw.pOp)), zeros(length(mrw.pOp)) .+ mrw.cov));
			# Update parameter values:
			for pI in 1:length(mrw.pOp)
				r0[pI] = p[mrw.pOp[pI]]; # Save previous value
				p[mrw.pOp[pI]] *= (mrw.M .^ rI[pI]); # Update value
				# Exclude values outside regime of exploration:
				if p[mrw.pOp[pI]] < (10.0 ^ mrw.pMin[pI])
					p[mrw.pOp[pI]] = (10.0 ^ mrw.pMin[pI])
				elseif p[mrw.pOp[pI]] > (10.0 ^ mrw.pMax[pI])
					p[mrw.pOp[pI]] = (10.0 ^ mrw.pMax[pI])
				end
			end
			# Calculate new simulated data + MSE:
			mse1 = fn.MSE(mySIM(fn,mm,p,d),d.Xe);
			# Evaluate if accept new parameter values or not:
			if(rand() < exp((mse0 - mse1) / mrwT[i]) || isnan(mse0))
				# If yes, update "reference" system
				mse0 = mse1;
			else
				# If not, revert to previous parameter values
				for pI in 1:length(mrw.pOp)
					p[mrw.pOp[pI]] = r0[pI];
				end
			end
			# Print results (either every iteration, prtW==1, or in the last iteration):
			if(mrw.prtW==1 || i==mrw.iter)
				writedlm(io, [vcat(ruN,i,mse0,[p[i] for i in mrw.pOp])],'\t')
			end
		end
		println("I: MSE = ",mse0)
		# TO DO: Print output simulated data
	end
end
