function [posPV, atLimit, changed] = q_limit_switch(V, Qg, posPV, atLimit, QgMin, QgMax, Vset, tol)
%Q_LIMIT_SWITCH  PV <-> PQ bus-type switching for generator reactive limits.
%
%   [posPV, atLimit, changed] = q_limit_switch(V, Qg, posPV, atLimit, QgMin, QgMax, Vset, tol)
%
%   V, Qg        : N x 1 |V| (pu) and generator reactive output (pu)
%   posPV        : N x 1 logical, buses currently regulating voltage
%   atLimit      : N x 1, +1 held at QgMax, -1 held at QgMin, 0 otherwise
%   QgMin, QgMax : N x 1 limits (pu); use -Inf/+Inf for unlimited buses
%   Vset         : N x 1 voltage set-points (pu)
%   tol          : struct with fields Q (pu) and V (pu)
%
%   Rules (standard power-flow practice):
%     PV bus,  Qg > QgMax            -> PQ at QgMax   (atLimit = +1)
%     PV bus,  Qg < QgMin            -> PQ at QgMin   (atLimit = -1)
%     at QgMax, |V| > Vset           -> back to PV    (it could lower Q)
%     at QgMin, |V| < Vset           -> back to PV    (it could raise Q)
%   changed marks every bus whose type changed.

toMax  = posPV & Qg > QgMax + tol.Q;
toMin  = posPV & Qg < QgMin - tol.Q;
backHi = atLimit == +1 & V > Vset + tol.V;
backLo = atLimit == -1 & V < Vset - tol.V;

posPV(toMax | toMin) = false;
atLimit(toMax) = +1;
atLimit(toMin) = -1;
back = backHi | backLo;
posPV(back) = true;
atLimit(back) = 0;
changed = toMax | toMin | back;
end
