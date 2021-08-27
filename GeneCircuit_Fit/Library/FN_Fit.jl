## ODE FITTING PIPELINE - Gene Regulatory Circuit models
#  ODE simulations & other functions
#	Mariana Gómez-Schiavon
#	August, 2021
#		Julia v.1.5.3
#		Local Environment: "./TOOLS/GeneCircuit_Fit/Project.toml"

module fn
	# Required libraries
	using DifferentialEquations
	using Distributions

	# Steady state function for a given system
	# INPUT: syst - Handle for the ODE system (@ode_def)
	#        p    - Dictionary function with the ODE parameters & values
	#        x0   - Vector of initial state of the ODE system
	#        rtol - Tolerance value for ODE solver
	# OUTPUT:     - Vector of steady state of the ODE system
	function SS(syst, p, x0, rtol)
		pV = [p[eval(Meta.parse(string(":",i)))] for i in syst.sys.ps];
		tS = 0;
		dXrm = 1;
		while(dXrm > rtol)
			ss = try
				solve(ODEProblem(syst,x0,1e6,pV); reltol=rtol,save_everystep = false);
			catch
				try
					solve(ODEProblem(syst,x0,1e6,pV),alg_hints=[:stiff]; reltol=rtol,save_everystep = false);
				catch err
					println("WARNING: Error in ODE simulation: <<",err,">>. ss --> NaN")
					x0 = zeros(length(syst.syms)).+NaN;
					break
				end
			end;
			dXrm = maximum(abs.(big.(ss(1e6))-big.(ss(1e6-0.01)))./big.(ss(1e6)));
			x0 = ss(1e6);
			tS += 1e6;
			if(tS>=1e12)
				println("WARNING: Maximum iteration reached (simulated time 1e12). Max relative Delta: ",dXrm)
				break
			end
		end
		return x0
	end;

	# ODE dynamics for a given system
	# INPUT: syst - Handle for the ODE system (@ode_def)
	#        p    - Dictionary function with the ODE parameters & values
	#        x0   - Vector of initial state of the ODE system
	#        tspan- Time to simulate
	#        rtol - Tolerance value for ODE solver
	# OUPUT: xD   - ODE system solver
	function Dyn(syst, p, x0, tspan, rtol)
		pV = [p[eval(Meta.parse(string(":",i)))] for i in syst.sys.ps];
		xD = try
			solve(ODEProblem(syst,x0,tspan,pV); reltol=rtol);
		catch
			try
				solve(ODEProblem(syst,x0,tspan,pV),alg_hints=[:stiff]; reltol=rtol);
			catch err
				println("WARNING: Error in ODE simulation: <<",err,">>. ss --> NaN")
				zeros(length(syst.syms)).+NaN;
			end
		end;
		return xD
	end;
end
