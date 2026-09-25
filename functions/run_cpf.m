function res = run_cpf(BusData, BranchData, opts)
%RUN_CPF  Continuation power flow (predictor-corrector, local parametrization).
%
%   res = run_cpf(BusData, BranchData, MonitorBus=14, Pfixed=..., ...)
%
%   Traces the P-V (nose) curve of the load/generation pattern in BusData:
%
%       P_inj(lambda) = Pfixed + lambda * (PG - PL)       (all buses, pu)
%       Q_inj(lambda) = Qfixed + lambda * (QG - QL)       (PQ buses)
%
%   lambda = 1 is the base case; lambda = 0 is no load. Pfixed/Qfixed model
%   DER (PV / wind farms) at constant output, i.e. NOT scaled with load.
%
%   Three phases (Crow, "Computational Methods for Electric Power Systems"):
%     1. lambda as continuation parameter, step StepLoad, until the
%        corrector fails near the nose;
%     2. |V| of MonitorBus as continuation parameter, step StepVoltage,
%        round the nose until lambda <= SwitchBackFraction * lambda_switch
%        (or, if |V| stops being a usable parameter past the nose, which
%        happens on radial feeders, until the corrector stalls);
%     3. lambda again (decreasing) down the lower branch until lambda < 0
%        or the corrector fails.
%
%   The state vector is z = [delta(non-slack); |V|(PQ); lambda], and the
%   corrector solves the augmented system  [f(z); z_p - z_p_pred] = 0  with
%   a full Newton loop.
%
%   Generator reactive limits (EnforceQLimits = true)
%     Limits are BusData columns 13 (Qmax) and 14 (Qmin), pu; a PV bus with
%     both equal to 0 is treated as unlimited. After every accepted step
%     q_limit_switch checks PV -> PQ (Qg outside limits) and PQ -> PV
%     (|V| back across the set-point). A crossing is located by halving the
%     step until it is <= EventTolerance, then the bus type is switched and
%     the point re-corrected with the same continuation parameter. The
%     starting point is solved with nrlf_qlim. Slack Q is not limited.
%     A limit hit can end the curve abruptly (limit-induced bifurcation):
%     the lambda-parametrized corrector then fails and phase 2 takes over.
%
%   Name-value options
%     MonitorBus          PQ bus used as continuation parameter in phase 2
%                         (must stay PQ; it cannot be a generator bus)
%     Pfixed, Qfixed      N x 1 constant injections (pu), default zeros
%     StepLoad            lambda step in phases 1 and 3        (0.1)
%     StepVoltage         |V| step in phase 2 (pu)             (0.005)
%     SwitchBackFraction  phase 2 -> 3 switch criterion        (0.75)
%     Tolerance           corrector mismatch tolerance (pu)    (1e-8)
%     MaxIterations       corrector Newton iterations          (25)
%     MaxSteps            safety cap on continuation steps     (5000)
%     MinStep             smallest step after halving (phase 2)(1e-5)
%     EnforceQLimits      generator Q limits on/off            (false)
%     EventTolerance      step length locating a Q-limit hit   (1e-4)
%
%   Output struct res
%     lambda (M x 1), V (N x M, pu), delta (N x M, rad), phase (M x 1),
%     Qg (N x M, generator Q in pu; NaN at non-generator buses),
%     atLimit (N x M, +1 at Qmax, -1 at Qmin, 0 regulating),
%     events (table: Point, Lambda, Bus, Action),
%     lambdaMax, VmonAtMax, idxMax, lambdaSwitch, monitorBus, Ybus,
%     posSL, posPV (initial types), qLimits (struct QgMin/QgMax/Vset)
%
%   Assumptions / limitations: balanced positive-sequence RMS model,
%   constant-power loads, slack bus absorbs all mismatch.

arguments
    BusData double
    BranchData double
    opts.MonitorBus (1,1) double {mustBeInteger, mustBePositive}
    opts.Pfixed double = []
    opts.Qfixed double = []
    opts.StepLoad (1,1) double {mustBePositive} = 0.1
    opts.StepVoltage (1,1) double {mustBePositive} = 0.005
    opts.SwitchBackFraction (1,1) double {mustBePositive} = 0.75
    opts.Tolerance (1,1) double {mustBePositive} = 1e-8
    opts.MaxIterations (1,1) double {mustBeInteger, mustBePositive} = 25
    opts.MaxSteps (1,1) double {mustBeInteger, mustBePositive} = 5000
    opts.MinStep (1,1) double {mustBePositive} = 1e-5
    opts.EnforceQLimits (1,1) logical = false
    opts.EventTolerance (1,1) double {mustBePositive} = 1e-4
end

N = size(BusData,1);
Pfixed = opts.Pfixed; if isempty(Pfixed), Pfixed = zeros(N,1); end
Qfixed = opts.Qfixed; if isempty(Qfixed), Qfixed = zeros(N,1); end
Pfixed = Pfixed(:); Qfixed = Qfixed(:);

Ybus  = Calculate_Ybus(BusData, BranchData);
posSL = BusData(:,4) == 3;
posPV0 = BusData(:,4) == 2;
isGen = posPV0 | posSL;
if ~(opts.MonitorBus <= N && BusData(opts.MonitorBus,4) == 0)
    error('run_cpf:MonitorBusNotPQ', ...
        'MonitorBus %d must be a PQ bus (type 0); its |V| is the phase-2 continuation parameter.', opts.MonitorBus);
end

Psch = BusData(:,9)  - BusData(:,7);    % load/generation direction vector
Qsch = BusData(:,10) - BusData(:,8);
Qsch(isGen) = -BusData(isGen,8);        % generator Q is a result, not a schedule

% Generator reactive limits (pu); +-Inf = unlimited
lim.Vset  = BusData(:,5);
lim.QgMax = inf(N,1);
lim.QgMin = -inf(N,1);
if opts.EnforceQLimits
    limited = posPV0 & ~(BusData(:,13) == 0 & BusData(:,14) == 0);
    lim.QgMax(limited) = BusData(limited,13);
    lim.QgMin(limited) = BusData(limited,14);
end
qTol = struct('Q', 1e-8, 'V', 1e-8);
maxEvents = 100;

% Near the nose the augmented matrix is (by design) close to singular in
% the lambda parametrization; that is detected via non-convergence.
w = warning('off', 'MATLAB:nearlySingularMatrix');
w(2) = warning('off', 'MATLAB:singularMatrix');
cleanup = onCleanup(@() warning(w)); %#ok<NASGU>

% ---- starting point: load flow at lambda = 0 (DER only) -----------------
V = ones(N,1);
d = zeros(N,1);
V(posSL|posPV0) = lim.Vset(posSL|posPV0);
[V, d, ok, posPV, atLimit] = nrlf_qlim(V, d, Ybus, Pfixed, Qfixed, posSL, posPV0, zeros(N,1), lim, ...
    Tolerance=opts.Tolerance, MaxIterations=opts.MaxIterations);
if ~ok
    error('run_cpf:NoStartingPoint', 'Load flow at lambda = 0 did not converge.');
end

% Bus-type dependent index sets (rebuilt after every Q-limit switch)
isAng = []; isMag = []; nA = 0; K = []; idxLambda = 0; idxVmon = 0; QgFix = [];
rebuild();
z = pack(0);

cap = opts.MaxSteps + 2*maxEvents + 1;
lamHist = nan(cap,1); phaseHist = nan(cap,1);
VHist = nan(N,cap); dHist = nan(N,cap); QgHist = nan(N,cap); limHist = zeros(N,cap);
evPoint = zeros(0,1); evLambda = zeros(0,1); evBus = zeros(0,1); evAction = strings(0,1);
nPts = 0;
record(z, 1);

% ---- continuation --------------------------------------------------------
phase = 1; s = +1; sigma = opts.StepLoad;
lambdaSwitch = NaN;
for step = 1:opts.MaxSteps
    p = paramIndex(phase);
    t = tangent(z, p, s);
    ok = all(isfinite(t));
    if ok
        [zNew, ok] = corrector(z + sigma*t, p);
    end

    % --- Q-limit event: locate, switch bus types, re-correct -------------
    if ok && opts.EnforceQLimits && ~(phase == 3 && zNew(end) < 0)
        [~, ~, changed] = limitCheck(zNew);
        if any(changed)
            if sigma > opts.EventTolerance
                sigma = sigma/2;
                continue
            end
            if numel(evBus) >= maxEvents
                warning('run_cpf:TooManyEvents', 'More than %d Q-limit events; curve truncated.', maxEvents);
                break
            end
            z = zNew; record(z, phase);
            [okSwitch, z] = switchTypes(z, phase);
            if ~okSwitch
                warning('run_cpf:SwitchFailed', ...
                    'Re-correction after a Q-limit switch failed at lambda = %.4f; curve truncated.', z(end));
                break
            end
            record(z, phase);
            sigma = baseStep(phase);
            continue
        end
    end

    switch phase
        case 1
            if ~ok
                lambdaSwitch = z(end);
                phase = 2; s = -1; sigma = opts.StepVoltage;
                continue
            end
            z = zNew; record(z, 1);
        case 2
            if ~ok
                sigma = sigma/2;
                if sigma >= opts.MinStep
                    continue
                end
                if z(end) < max(lamHist(1:nPts)) && any(phaseHist(1:nPts) == 2)
                    % Already past the nose, but |V_mon| has stopped being a
                    % good parameter (it can fold on the lower branch of a
                    % radial feeder): continue with lambda instead.
                    phase = 3; s = -1; sigma = opts.StepLoad;
                    continue
                end
                warning('run_cpf:Phase2Stalled', ...
                    'Voltage-parametrized corrector failed at lambda = %.4f; curve truncated.', z(end));
                break
            end
            z = zNew; record(z, 2);
            if z(end) <= opts.SwitchBackFraction*lambdaSwitch
                phase = 3; s = -1; sigma = opts.StepLoad;
            end
        case 3
            if ~ok || zNew(end) < 0
                break
            end
            z = zNew; record(z, 3);
    end
end
if step == opts.MaxSteps
    warning('run_cpf:MaxSteps', 'MaxSteps (%d) reached before the curve was completed.', opts.MaxSteps);
end

res.lambda  = lamHist(1:nPts);
res.V       = VHist(:,1:nPts);
res.delta   = dHist(:,1:nPts);
res.phase   = phaseHist(1:nPts);
res.Qg      = QgHist(:,1:nPts);
res.atLimit = limHist(:,1:nPts);
res.events  = table(evPoint, evLambda, evBus, evAction, ...
    'VariableNames', {'Point', 'Lambda', 'Bus', 'Action'});
[res.lambdaMax, res.idxMax] = max(res.lambda);
res.VmonAtMax    = res.V(opts.MonitorBus, res.idxMax);
res.lambdaSwitch = lambdaSwitch;
res.monitorBus   = opts.MonitorBus;
res.Ybus         = Ybus;
res.posSL = posSL; res.posPV = posPV0;
res.qLimits = lim;
res.qLimitsEnforced = opts.EnforceQLimits;

% ======================================================================
    function rebuild()
        isAng = ~posSL;
        isMag = ~(posSL | posPV);
        nA = nnz(isAng);
        magBuses = find(isMag);
        K = [Psch(isAng); Qsch(isMag)];
        idxLambda = nA + numel(magBuses) + 1;
        idxVmon = nA + find(magBuses == opts.MonitorBus, 1);
        QgFix = zeros(N,1);                         % generator Q held at a limit
        QgFix(atLimit == +1) = lim.QgMax(atLimit == +1);
        QgFix(atLimit == -1) = lim.QgMin(atLimit == -1);
    end

    function p = paramIndex(ph)
        if ph == 2, p = idxVmon; else, p = idxLambda; end
    end

    function st = baseStep(ph)
        if ph == 2, st = opts.StepVoltage; else, st = opts.StepLoad; end
    end

    function zz = pack(lam)
        zz = [d(isAng); V(isMag); lam];
    end

    function [Vf, df, lam] = unpack(zz)
        Vf = V; df = d;
        df(isAng) = zz(1:nA);
        Vf(isMag) = zz(nA+1:end-1);
        lam = zz(end);
    end

    function [F, A, Qc] = augmented(zz, pIdx)
        % F: power mismatch, A: [J_|V| -K; e_p], Qc: calculated Q injection
        [Vf, df, lam] = unpack(zz);
        [~, Qc, ~, ~, ~, ~, F] = Calculate_PcalcQcalc(Vf, df, Ybus, ...
            Pfixed + lam*Psch, Qfixed + QgFix + lam*Qsch, posSL, posPV);
        if nargout > 1
            [~, J] = Jacobian_NRLF(Ybus, Vf, df, posSL, posPV);
            J(:, nA+1:end) = J(:, nA+1:end) ./ transpose(Vf(isMag));  % d/d|V|
            ep = zeros(1, numel(zz)); ep(pIdx) = 1;
            A = [J, -K; ep];
        end
    end

    function t = tangent(zz, pIdx, sgn)
        [~, A] = augmented(zz, pIdx);
        rhs = zeros(numel(zz),1); rhs(end) = sgn;
        t = A\rhs;
    end

    function [zz, ok] = corrector(zz, pIdx)
        target = zz(pIdx);
        ok = false;
        for it = 1:opts.MaxIterations + 1
            [F, A] = augmented(zz, pIdx);
            if ~all(isfinite(F)), return; end
            if max(abs(F)) <= opts.Tolerance && abs(zz(pIdx)-target) <= opts.Tolerance
                ok = all(zz(nA+1:end-1) > 0);
                return
            end
            if it > opts.MaxIterations, return; end
            dz = A \ [F; target - zz(pIdx)];
            if ~all(isfinite(dz)), return; end
            zz = zz + dz;
        end
    end

    function Qg = generatorQ(zz)
        [~, ~, Qc] = augmented(zz, idxLambda);
        lam = zz(end);
        Qg = nan(N,1);
        Qg(isGen) = Qc(isGen) - Qfixed(isGen) - lam*Qsch(isGen);
    end

    function [newPV, newLim, changed] = limitCheck(zz)
        [Vf, ~, ~] = unpack(zz);
        [newPV, newLim, changed] = q_limit_switch(Vf, generatorQ(zz), posPV, atLimit, ...
            lim.QgMin, lim.QgMax, lim.Vset, qTol);
    end

    function [ok, zz] = switchTypes(zz, ph)
        % Apply the switches found at zz and re-solve at the same value of
        % the continuation parameter. Repeats if the re-solved point
        % triggers further switches (simultaneous limits).
        for round = 1:10
            [newPV, newLim, changed] = limitCheck(zz);
            if ~any(changed)
                ok = true;
                return
            end
            [Vf, df, lam] = unpack(zz);
            for b = find(changed).'
                evPoint(end+1,1) = nPts; evLambda(end+1,1) = lam; evBus(end+1,1) = b; %#ok<AGROW>
                if newPV(b)
                    evAction(end+1,1) = "PQ->PV";           %#ok<AGROW>
                    Vf(b) = lim.Vset(b);
                elseif newLim(b) > 0
                    evAction(end+1,1) = "PV->PQ at Qmax";   %#ok<AGROW>
                else
                    evAction(end+1,1) = "PV->PQ at Qmin";   %#ok<AGROW>
                end
            end
            V = Vf; d = df;
            posPV = newPV; atLimit = newLim;
            rebuild();
            zz = pack(lam);
            if ph == 2
                % keep |V_mon| at the value it had before the switch
                [zz, ok] = corrector(zz, idxVmon);
            else
                [zz, ok] = corrector(zz, idxLambda);
            end
            if ~ok, return; end
        end
        ok = false;
    end

    function record(zz, ph)
        [Vf, df, lam] = unpack(zz);
        V = Vf; d = df;             % keep PV/slack entries, warm start
        nPts = nPts + 1;
        lamHist(nPts) = lam; phaseHist(nPts) = ph;
        VHist(:,nPts) = Vf;  dHist(:,nPts) = df;
        QgHist(:,nPts) = generatorQ(zz);
        limHist(:,nPts) = atLimit;
    end
end
