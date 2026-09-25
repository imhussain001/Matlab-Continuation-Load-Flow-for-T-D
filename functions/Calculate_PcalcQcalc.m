function [Pcalc,Qcalc,dP,dQ,dPred,dQred,dPdQred] = Calculate_PcalcQcalc(V,d,Ybus,Psch,Qsch,posSL,posPV)
%CALCULATE_PCALCQCALC  Injected P/Q and mismatches (pu).
%
%   V, d       : N x 1 voltage magnitude (pu) and angle (rad)
%   Psch, Qsch : N x 1 scheduled net injections (generation - load, pu)
%   posSL/posPV: N x 1 logical masks of slack / PV buses
%
%   dPdQred = [dP at non-slack buses; dQ at PQ buses] is the NR mismatch.

Vc = V.*exp(1i*d);
Scalc = conj(Vc).*(Ybus*Vc);   % = conj(S_injected)
Pcalc = real(Scalc);
Qcalc = -imag(Scalc);
dP = Psch - Pcalc;
dQ = Qsch - Qcalc;
dPred = dP(posSL==0);
dQred = dQ(posSL+posPV==0);
dPdQred = [dPred; dQred];
end
