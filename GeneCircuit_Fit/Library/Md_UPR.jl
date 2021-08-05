## ODE FITTING PIPELINE - Gene Regulatory Circuit models
# ODE model - Unfolded protein response as modeled in Pincus et al. (2010)
#	Mariana Gómez-Schiavon
#	August, 2021
#		Julia v.1.5.3
#		Local Environment: "./TOOLS/GeneCircuit_Fit/Project.toml"

module mm
	# Required libraries
	using DifferentialEquations
	using ParameterizedFunctions

	# ODE system
	myODE = @ode_def begin
		dU   = sU - (cB * U * B) + (gUB * UB) - (cD * U) + (cE * E * Ud) + (gB * UB)
		dUB  = (cB * U * B) - (gUB * UB) - (cD * UB) + (cE * E * UdB) - (gB * UB) - (gF * UB)
		dUd  = - (cB * Ud * B) + (gUB * UdB) + (cD * U) - (cE * E * Ud) + (gB * UdB)
		dUdB = (cB * Ud * B) - (gUB * UdB) + (cD * UB) - (cE * E * UdB) - (gB * UdB)
		dI   = mI - (gI * I) - (cBI * B * I) - (cA * (U + Ud + UdB) * I) + (gIB * IB) + (gIA * ((kI ^ nI)/((kI ^ nI) + (Ia ^ nI))) * Ia)
		dIB  = - (gI * IB) + (cBI * B * I) - (gIB * IB)
		dIa  = - (gI * Ia) + (cA * (U + Ud + UdB) * I) - (gIA * ((kI ^ nI)/((kI ^ nI) + (Ia ^ nI))) * Ia)
		dHu  = bHu - (gHs * Hu) - ((bHo * min(I,Hu)) + (bHs * min(Ia,Hu)))
		dHs  = - (gHs * Hs) + ((bHo * min(I,Hu)) + (bHs * min(Ia,Hu)))
		dBm  = (bBm * (1 + (nB * (Hs ^ 2)/(a0 + (a1 * Hs) + (Hs ^ 2))))) - (gBm * Bm)
		dB   = (bB * Bm) - (gB * B) - (cB * U * B) + (gUB * UB) - (cB * Ud * B) + (gUB * UdB) + (gF * UB) - (cBI * B * I) + (gIB * IB) + (gI * IB)
		dEm  = (bEm * (1 + (nE * (Hs ^ 2)/(a0 + (a1 * Hs) + (Hs ^ 2))))) - (gEm * Em)
		dE   = (bE * Em) - (gE * E)
		dRm  = - (gHs * Rm) + ((bHo * min(I,Hu)) + (bHs * min(Ia,Hu)))
		dRs  = (0.008333 * Rm) - (0.00003472 * Rs)
	end sU cB cBI gUB cD cE gB gF mI gI cA gIB gIA kI nI bHu gHs bHs bHo bBm nB a0 a1 gBm bB gB bEm nE gEm bE gE uT;
end
