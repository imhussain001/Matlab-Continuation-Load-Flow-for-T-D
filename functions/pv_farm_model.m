function pv = pv_farm_model(opts)
%PV_FARM_MODEL  Steady-state PV plant: single-diode panel -> array -> boost -> VSI.
%
%   pv = pv_farm_model(PTarget_W=1.5e6, PowerFactor=0.8071, ...)
%
%   Panel: single-diode model with ideality factor a = 2 (as in the
%   original pvmpdel3.m), solved point-by-point with Newton-Raphson.
%   Array: N_series panels per string to reach VArrayTarget_V at MPP,
%   N_parallel strings to reach PTarget_W. The plant runs at MPP.
%   Reactive power: fixed power factor, Q = P*tan(acos(pf)), injected.
%
%   Output fields (SI units)
%     P_W, Q_W                 plant active / reactive output
%     N_series, N_parallel, N_panels
%     Vmp_V, Imp_A, Pmp_W      single-panel MPP
%     VArray_V, IArray_A       array MPP
%     DutyCycle, Vdc_V, ModIndex   boost converter / inverter operating point
%     curve.V, curve.I, curve.P    single-panel I-V and P-V curves
%
%   Legacy notes: the first pvmpdel3.m version modelled ONE string
%   (P = 3.1 kW -> 3.1e-5 pu on 100 MVA), so the "PV" cases had no
%   measurable effect. PVData.mat and pv_data.mat also disagreed (Q =
%   1.50 kvar vs 1.92 kvar). The pf 0.8071 comment "cos(45 deg)" was wrong:
%   cos(45 deg) = 0.7071; 0.8071 gives Q/P = 0.73.

arguments
    opts.PTarget_W (1,1) double {mustBePositive} = 1.5e6
    opts.PowerFactor (1,1) double {mustBeInRange(opts.PowerFactor,0,1,"exclude-lower")} = 0.8071
    opts.Irradiance_Wm2 (1,1) double {mustBeNonnegative} = 1000
    opts.CellTemp_C (1,1) double = 25
    opts.VArrayTarget_V (1,1) double {mustBePositive} = 440
    opts.VdcTarget_V (1,1) double {mustBePositive} = 650
    opts.VGrid_V (1,1) double {mustBePositive} = 400
end

% Panel datasheet (72-cell module)
kB   = 1.38065e-23;   % Boltzmann constant (J/K)
qe   = 1.602e-19;     % electron charge (C)
Iscn = 8.21;          % short-circuit current at STC (A)
Vocn = 32.9;          % open-circuit voltage at STC (V)
Ns   = 72;            % cells in series per panel
Gn   = 1000;          % STC irradiance (W/m^2)
Rs   = 0.15;          % series resistance (ohm)
Rp   = 500;           % parallel resistance (ohm)
a    = 2;             % diode ideality factor

T   = opts.CellTemp_C + 273.15;
Vt  = Ns*kB*T/qe;                           % panel thermal voltage (V)
I0  = Iscn/(exp(Vocn/(a*Vt)) - 1);          % saturation current (A)
Iph = (opts.Irradiance_Wm2/Gn)*Iscn;        % photocurrent (A)

V = linspace(0, Vocn, 400);
I = zeros(size(V));
for k = 1:numel(V)
    Ik = Iph;
    for it = 1:50
        e  = exp((V(k) + Ik*Rs)/(a*Vt));
        f  = Ik - Iph + I0*(e - 1) + (V(k) + Ik*Rs)/Rp;
        df = 1 + I0*Rs/(a*Vt)*e + Rs/Rp;
        dI = f/df;
        Ik = Ik - dI;
        if abs(dI) < 1e-12, break; end
    end
    I(k) = Ik;
end
I = max(I, 0);
P = V.*I;
[Pmp, idx] = max(P);

pv.Vmp_V = V(idx);
pv.Imp_A = I(idx);
pv.Pmp_W = Pmp;
pv.curve = struct('V', V, 'I', I, 'P', P);

% Array sizing
pv.N_series   = ceil(opts.VArrayTarget_V/pv.Vmp_V);
pv.N_parallel = ceil(ceil(opts.PTarget_W/Pmp)/pv.N_series);
pv.N_panels   = pv.N_series*pv.N_parallel;
pv.VArray_V   = pv.N_series*pv.Vmp_V;
pv.IArray_A   = pv.N_parallel*pv.Imp_A;
pv.P_W        = pv.VArray_V*pv.IArray_A;
pv.Q_W        = pv.P_W*tan(acos(opts.PowerFactor));
pv.PowerFactor = opts.PowerFactor;

% Boost converter and inverter operating point (informational only; not
% used by the load flow)
pv.DutyCycle = min(max(1 - pv.VArray_V/opts.VdcTarget_V, 0.3), 0.6);
pv.Vdc_V     = pv.VArray_V/(1 - pv.DutyCycle);
pv.ModIndex  = min(max(opts.VGrid_V*sqrt(2)/pv.Vdc_V, 0.85), 0.95);
end
