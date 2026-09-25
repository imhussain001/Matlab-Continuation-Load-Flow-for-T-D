function [Sder_pu, info] = der_injection(scenario, cfg, N, Sbase_MVA)
%DER_INJECTION  Constant complex DER injection vector for a scenario.
%
%   [Sder_pu, info] = der_injection(scenario, cfg, N, Sbase_MVA)
%
%   scenario : 'Base' | 'PV' | 'Wind' | 'PV+Wind'
%   Returns an N x 1 complex vector (pu, generation positive) with the
%   selected plants at cfg.derBus, and info with the plant models and
%   totals in MW / Mvar.

hasPV   = any(strcmpi(scenario, {'PV', 'PV+Wind'}));
hasWind = any(strcmpi(scenario, {'Wind', 'PV+Wind'}));
if ~hasPV && ~hasWind && ~strcmpi(scenario, 'Base')
    error('der_injection:UnknownScenario', 'Unknown scenario "%s".', scenario);
end

P_W = 0; Q_W = 0;
info.pv = []; info.wind = [];
if hasPV
    info.pv = pv_farm_model(cfg.pv{:});
    P_W = P_W + info.pv.P_W;  Q_W = Q_W + info.pv.Q_W;
end
if hasWind
    info.wind = wind_farm_model(cfg.wind{:});
    P_W = P_W + info.wind.P_W;  Q_W = Q_W + info.wind.Q_W;
end

Sder_pu = complex(zeros(N,1));
Sder_pu(cfg.derBus) = complex(P_W, Q_W)/(Sbase_MVA*1e6);
info.P_MW = P_W/1e6;
info.Q_Mvar = Q_W/1e6;
end
