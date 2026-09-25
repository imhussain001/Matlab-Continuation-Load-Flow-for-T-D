function [V,d,Pcalc,Qcalc,dP,dQ,dPred,dQred,dPdQred,it,converged] = NRLF(V,d,Ybus,Psch,Qsch,posSL,posPV,Toler,Max_It)
%NRLF  Newton-Raphson load flow (polar form).
%
%   [V,d,...,it,converged] = NRLF(V,d,Ybus,Psch,Qsch,posSL,posPV,Toler,Max_It)
%
%   V, d       : initial guess, N x 1 (pu, rad)
%   Psch, Qsch : scheduled net injections (pu)
%   Toler      : max-abs mismatch tolerance (pu)
%   Max_It     : maximum Newton iterations
%
%   converged is true only if the final mismatch is finite and <= Toler.
%   (Previously a NaN mismatch was silently treated as converged, and
%   callers had to infer success from the iteration count.)
%   PV-bus reactive limits are NOT enforced.

it = 0;
N = length(V);
while true
    [Pcalc,Qcalc,dP,dQ,dPred,dQred,dPdQred] = Calculate_PcalcQcalc(V,d,Ybus,Psch,Qsch,posSL,posPV);
    converged = all(isfinite(dPdQred)) && max(abs(dPdQred)) <= Toler;
    if converged || it >= Max_It || ~all(isfinite(dPdQred))
        break;
    end
    it = it + 1;
    [~,J] = Jacobian_NRLF(Ybus,V,d,posSL,posPV);
    dddVred = J\dPdQred;              % [d(delta); d|V|/|V|]
    dddV = 1 - [posSL; posSL+posPV];
    dddV(dddV==1) = dddVred;
    d = d + dddV(1:N);
    V = V.*(1 + dddV(N+1:2*N));
end
end
