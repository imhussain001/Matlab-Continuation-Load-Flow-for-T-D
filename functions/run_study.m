function out = run_study(systemName)
%RUN_STUDY  Solve the four scenarios (Base, PV, Wind, PV+Wind) of one system.
%
%   out = run_study('transmission')   or   run_study('distribution')
%
%   For every scenario:
%     1. PV / wind plant models -> constant P and Q injected at cfg.derBus
%     2. continuation power flow -> loadability (P-V) curve and lambda_max
%     3. load flow at base load (lambda = 1) -> voltages, losses and the
%        stability index (L-index for transmission, FVSI for distribution)
%
%   Returns
%     out.cfg, out.system       settings and network data
%     out.scenarios(k)          full results of scenario k
%     out.summary               one-row-per-scenario table (no NaN)
%     out.plants                PV plant and wind farm ratings
%
%   Nothing is plotted or saved here - run_all does that.

cfg = study_config(systemName);
isRadial = strcmp(cfg.system, 'ieee33');
if isRadial
    sys = ieee33_data();
    net = build_bibc_bcbv(sys.fromBus, sys.toBus, sys.z);
else
    sys = ieee14_data();
end
B = sys.BusData;
N = size(B,1);

nS = numel(cfg.scenarios);
sc = repmat(struct('name','','der',[],'cpf',[],'cpfNoQLim',[],'op',[],'index',[]), nS, 1);
for k = 1:nS
    name = cfg.scenarios{k};
    fprintf('  %s: solving %-8s ...\n', cfg.label, name);
    [Sder, der] = der_injection(name, cfg, N, sys.Sbase_MVA);

    cpf = run_cpf(B, sys.BranchData, 'MonitorBus', cfg.paramBus, ...
        'Pfixed', real(Sder), 'Qfixed', imag(Sder), cfg.cpf{:});
    op = operating_point(B, cpf, Sder, 1.0);
    if cpf.qLimitsEnforced && any(cpf.posPV)
        % Same case without generator Q limits, to show their effect
        sc(k).cpfNoQLim = run_cpf(B, sys.BranchData, 'MonitorBus', cfg.paramBus, ...
            'Pfixed', real(Sder), 'Qfixed', imag(Sder), cfg.cpf{:}, 'EnforceQLimits', false);
    end

    if isRadial
        % Independent backward/forward sweep at base load (also gives branch currents)
        [Vbfs, okBfs, ~, Ibr] = bfs_load_flow(net, sys.Sload - Sder);
        assert(okBfs, 'run_study:BfsFailed', 'BFS did not converge at base load (%s).', name);
        op.loss_MW = sum(abs(Ibr).^2 .* real(sys.z))*sys.Sbase_MVA;
        idx.name = 'FVSI';
        idx.value = fvsi_index(sys.fromBus, sys.toBus, sys.z, Vbfs, Ibr);
        idx.labels = sys.toBus;          % branch k ends at bus toBus(k)
    else
        % A generator held at its Q limit acts as a load bus in the L-index;
        % only the original load buses are reported so the scenarios line up.
        typeOp = B(:,4); typeOp(op.atLimit ~= 0) = 0;
        [Lv, loadBuses] = l_index(cpf.Ybus, op.Vc, typeOp);
        keep = B(loadBuses,4) == 0;
        idx.name = 'L-index';
        idx.value = Lv(keep);
        idx.labels = loadBuses(keep);
        op.loss_MW = real(sum(op.Vc .* conj(cpf.Ybus*op.Vc)))*sys.Sbase_MVA;
    end

    sc(k).name = name; sc(k).der = der; sc(k).cpf = cpf; sc(k).op = op; sc(k).index = idx;
end

out.cfg = cfg;
out.system = sys;
out.scenarios = sc;
out.summary = summary_table(out);
out.plants = plant_table(out);
end

% ------------------------------------------------------------------------
function op = operating_point(B, cpf, Sder, lambda)
% Load flow (same generator Q limits as the CPF) at the given loading,
% warm-started from the nearest upper-branch CPF point so that the stable
% high-voltage solution is found.
upper = find(cpf.phase == 1);
[~, j] = min(abs(cpf.lambda(upper) - lambda));
k = upper(j);
V = cpf.V(:,k); d = cpf.delta(:,k);
atLimit = cpf.atLimit(:,k);
posPV = cpf.posPV & atLimit == 0;
isGen = cpf.posPV | cpf.posSL;
Psch = real(Sder) + lambda*(B(:,9) - B(:,7));
Qother = imag(Sder) + lambda*(B(:,10) - B(:,8));
Qother(isGen) = imag(Sder(isGen)) - lambda*B(isGen,8);   % excludes generator Q
[V, d, ok, ~, atLimit, Qg] = nrlf_qlim(V, d, cpf.Ybus, Psch, Qother, cpf.posSL, posPV, atLimit, cpf.qLimits);
assert(ok, 'run_study:OperatingPoint', 'Load flow at lambda = %.2f did not converge.', lambda);
op.lambda = lambda;
op.V = V; op.delta = d; op.Vc = V.*exp(1i*d);
op.atLimit = atLimit; op.Qg = Qg;
end

function T = summary_table(out)
% One row per scenario, units in the column names, no NaN.
sc = out.scenarios; cfg = out.cfg; B = out.system.BusData;
baseLoad_MW = sum(B(:,7))*out.system.Sbase_MVA;
loadBus = B(:,4) == 0;
n = numel(sc);
T = table(string({sc.name}).', 'VariableNames', {'Scenario'});
T.PV_MW          = arrayfun(@(s) plant_value(s.der.pv,   'P_W'), sc)/1e6;
T.Wind_MW        = arrayfun(@(s) plant_value(s.der.wind, 'P_W'), sc)/1e6;
T.Renewable_Mvar = arrayfun(@(s) s.der.Q_Mvar, sc);
T.LambdaMax      = arrayfun(@(s) s.cpf.lambdaMax, sc);
if ~isempty(sc(1).cpfNoQLim)
    T.LambdaMax_noQlimits = arrayfun(@(s) s.cpfNoQLim.lambdaMax, sc);
end
T.LoadMargin_MW  = (T.LambdaMax - 1)*baseLoad_MW;
for b = cfg.plotBuses
    T.(sprintf('V%d_base_pu', b)) = arrayfun(@(s) s.op.V(b), sc);
end
Vmin = zeros(n,1); VminBus = zeros(n,1);
for k = 1:n
    v = sc(k).op.V; v(~loadBus) = inf;
    [Vmin(k), VminBus(k)] = min(v);
end
T.Vmin_base_pu = Vmin;
T.Vmin_bus = VminBus;
T.Losses_MW = arrayfun(@(s) s.op.loss_MW, sc);
T.([strrep(sc(1).index.name, '-', '') '_max']) = arrayfun(@(s) max(s.index.value), sc);
end

function T = plant_table(out)
% PV plant and wind farm ratings (identical in every scenario that uses them)
sc = out.scenarios; sys = out.system;
baseLoad_MW = sum(sys.BusData(:,7))*sys.Sbase_MVA;
pv = sc(find(arrayfun(@(s) ~isempty(s.der.pv), sc), 1)).der.pv;
wd = sc(find(arrayfun(@(s) ~isempty(s.der.wind), sc), 1)).der.wind;
Plant = ["PV plant"; "Wind farm"];
Size = [sprintf("%d panels (%d x %d)", pv.N_panels, pv.N_series, pv.N_parallel); ...
        sprintf("%d turbines x %.2f MW", wd.N_turbines, wd.PerTurbineP_W/1e6)];
P_MW = [pv.P_W; wd.P_W]/1e6;
Q_Mvar = [pv.Q_W; wd.Q_W]/1e6;
PowerFactor = P_MW./hypot(P_MW, Q_Mvar);
PercentOfBaseLoad = 100*P_MW/baseLoad_MW;
System = repmat(string(out.cfg.label), 2, 1);
Bus = repmat(out.cfg.derBus, 2, 1);
T = table(System, Plant, Bus, Size, P_MW, Q_Mvar, PowerFactor, PercentOfBaseLoad);
end

function v = plant_value(plant, field)
if isempty(plant), v = 0; else, v = plant.(field); end
end
