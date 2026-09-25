function [V, d, converged, posPV, atLimit, Qg] = nrlf_qlim(V, d, Ybus, Psch, Qother, posSL, posPV, atLimit, lim, opts)
%NRLF_QLIM  Newton-Raphson load flow with generator reactive-power limits.
%
%   [V,d,converged,posPV,atLimit,Qg] = nrlf_qlim(V,d,Ybus,Psch,Qother,posSL,posPV,atLimit,lim)
%
%   Qother  : N x 1 net reactive injection EXCLUDING the generator's own Q
%             (e.g. -QL + Q_DER), pu. At PQ buses without a generator this
%             is the whole scheduled Q.
%   posPV   : initial regulating buses; atLimit: initial limit state
%             (zeros(N,1) for a fresh solve, see q_limit_switch)
%   lim     : struct with QgMin, QgMax, Vset (N x 1, pu; +-Inf = unlimited)
%
%   Outer loop: solve NRLF, then switch PV <-> PQ with q_limit_switch
%   until no bus changes type. Qg is the generator reactive output (pu).
%
%   Name-value: Tolerance (1e-10), MaxIterations (30), MaxSwitchRounds (20)

arguments
    V double
    d double
    Ybus double
    Psch double
    Qother double
    posSL logical
    posPV logical
    atLimit double
    lim struct
    opts.Tolerance (1,1) double {mustBePositive} = 1e-10
    opts.MaxIterations (1,1) double {mustBeInteger, mustBePositive} = 30
    opts.MaxSwitchRounds (1,1) double {mustBeInteger, mustBePositive} = 20
end

tol = struct('Q', 1e-8, 'V', 1e-8);
converged = false;
Qg = nan(size(V));
for round = 1:opts.MaxSwitchRounds
    Qsch = Qother + (atLimit == +1).*finite0(lim.QgMax) + (atLimit == -1).*finite0(lim.QgMin);
    V(posPV) = lim.Vset(posPV);
    [V, d, ~, Qcalc, ~, ~, ~, ~, ~, ~, ok] = NRLF(V, d, Ybus, Psch, Qsch, posSL, posPV, ...
        opts.Tolerance, opts.MaxIterations);
    if ~ok
        return
    end
    Qg = Qcalc - Qother;
    [posPV, atLimit, changed] = q_limit_switch(V, Qg, posPV, atLimit, lim.QgMin, lim.QgMax, lim.Vset, tol);
    if ~any(changed)
        converged = true;
        return
    end
end
end

function x = finite0(x)
x(~isfinite(x)) = 0;
end
