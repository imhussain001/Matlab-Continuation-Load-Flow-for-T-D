function cfg = study_config(systemName)
%STUDY_CONFIG  All study parameters in one place - edit this file to change the study.
%
%   cfg = study_config('transmission')   % IEEE 14-bus  (alias 'ieee14')
%   cfg = study_config('distribution')   % IEEE 33-bus  (alias 'ieee33')
%
%   Every scenario (Base, PV, Wind, PV+Wind) of a system uses the same
%   renewable bus, plant ratings and CPF settings, so the scenarios differ
%   ONLY in which plants are connected.
%
%   Renewables are constant P/Q injections: they do NOT grow with the
%   load multiplier lambda.

switch lower(systemName)
    case {'transmission', 'ieee14'}
        cfg.system     = 'ieee14';
        cfg.label      = 'Transmission';
        cfg.title      = 'Transmission system (IEEE 14-bus)';
        cfg.derBus     = 14;                % PV + wind connection bus
        cfg.plotBuses  = 14;                % bus(es) shown in the loadability figure
        cfg.paramBus   = 14;                % PQ bus used by the CPF around the nose
        cfg.pv   = {'PTarget_W', 1.3e6, 'PowerFactor', 0.8071};
        cfg.wind = {'PTarget_W', 1.5e6, 'QTarget_W', 1.8e6, 'QOverP', 0.3};
        % Generator Q limits (BusData cols 13/14, generators at buses 2, 3, 6, 8)
        cfg.cpf  = {'StepLoad', 0.1, 'StepVoltage', 0.005, 'Tolerance', 1e-8, ...
                    'EnforceQLimits', true};
    case {'distribution', 'ieee33'}
        cfg.system     = 'ieee33';
        cfg.label      = 'Distribution';
        cfg.title      = 'Distribution system (IEEE 33-bus)';
        cfg.derBus     = 33;
        cfg.plotBuses  = [33 18];           % renewable bus and weakest bus
        % Bus 18 collapses first; |V33| folds back after the nose, so bus 18
        % is the better CPF parameter even though the plants are at bus 33.
        cfg.paramBus   = 18;
        cfg.pv   = {'PTarget_W', 1.5e6, 'PowerFactor', 0.8071};
        cfg.wind = {'PTarget_W', 2.5e6, 'QTarget_W', 2.2e6, 'PowerFactor', 0.85};
        % Loads are ~1e-3 pu on 100 MVA, hence the tight tolerance.
        % Only the substation is a source, so Q limits have no effect here.
        cfg.cpf  = {'StepLoad', 0.1, 'StepVoltage', 0.005, 'Tolerance', 1e-10, ...
                    'EnforceQLimits', true};
    otherwise
        error('study_config:UnknownSystem', ...
            'Unknown system "%s" (use ''transmission'' or ''distribution'').', systemName);
end

cfg.scenarios = {'Base', 'PV', 'Wind', 'PV+Wind'};
end
