function wf = wind_farm_model(opts)
%WIND_FARM_MODEL  Steady-state wind farm: aerodynamics -> generator -> boost -> VSI.
%
%   wf = wind_farm_model(PTarget_W=1.5e6, QTarget_W=1.8e6, QOverP=0.3)
%   wf = wind_farm_model(PTarget_W=2.5e6, QTarget_W=2.2e6, PowerFactor=0.85)
%
%   Per turbine: P_mech = 0.5 rho A Cp v^3, P_e = eta P_mech.
%   Reactive power per turbine is either QOverP*P_e or P_e*tan(acos(pf))
%   (give exactly one of QOverP / PowerFactor).
%   Farm size (legacy windss.m rule, kept deliberately):
%       N = max(ceil(PTarget/P_e), ceil(QTarget/Q_e))
%   so the turbine count can be driven by the reactive target. Set
%   QTarget_W = 0 to size on active power only.
%
%   Output fields (SI units)
%     P_W, Q_W, N_turbines, PerTurbineP_W, PerTurbineQ_W, PMech_W,
%     DutyCycle, Vdc_V, ModIndex

arguments
    opts.PTarget_W (1,1) double {mustBePositive} = 1.5e6
    opts.QTarget_W (1,1) double {mustBeNonnegative} = 0
    opts.QOverP double = []
    opts.PowerFactor double = []
    opts.WindSpeed_ms (1,1) double {mustBePositive} = 10
    opts.RotorRadius_m (1,1) double {mustBePositive} = 45
    opts.Cp (1,1) double {mustBeInRange(opts.Cp,0,0.593)} = 0.4    % Betz limit
    opts.AirDensity_kgm3 (1,1) double {mustBePositive} = 1.225
    opts.Efficiency (1,1) double {mustBeInRange(opts.Efficiency,0,1)} = 0.9
    opts.VGen_V (1,1) double {mustBePositive} = 400
    opts.VdcTarget_V (1,1) double {mustBePositive} = 650
    opts.VGrid_V (1,1) double {mustBePositive} = 400
end

if isempty(opts.QOverP) == isempty(opts.PowerFactor)
    error('wind_farm_model:QSpec', 'Specify exactly one of QOverP or PowerFactor.');
end

A = pi*opts.RotorRadius_m^2;
wf.PMech_W = 0.5*opts.AirDensity_kgm3*A*opts.Cp*opts.WindSpeed_ms^3;
wf.PerTurbineP_W = opts.Efficiency*wf.PMech_W;
if ~isempty(opts.QOverP)
    wf.PerTurbineQ_W = opts.QOverP*wf.PerTurbineP_W;
else
    wf.PerTurbineQ_W = wf.PerTurbineP_W*tan(acos(opts.PowerFactor));
end

nP = ceil(opts.PTarget_W/wf.PerTurbineP_W);
nQ = 0;
if opts.QTarget_W > 0 && wf.PerTurbineQ_W > 0
    nQ = ceil(opts.QTarget_W/wf.PerTurbineQ_W);
end
wf.N_turbines = max(nP, nQ);
wf.P_W = wf.N_turbines*wf.PerTurbineP_W;
wf.Q_W = wf.N_turbines*wf.PerTurbineQ_W;

% Converter operating point (informational only)
wf.DutyCycle = min(max(1 - opts.VGen_V/opts.VdcTarget_V, 0.3), 0.6);
wf.Vdc_V     = opts.VGen_V/(1 - wf.DutyCycle);
wf.ModIndex  = min(max(opts.VGrid_V*sqrt(2)/wf.Vdc_V, 0.85), 0.95);
end
