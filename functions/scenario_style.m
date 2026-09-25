function [col, mk] = scenario_style(scenario)
%SCENARIO_STYLE  Colour and marker of a scenario - identical in every figure.
switch scenario
    case 'Base',    col = [0.15 0.15 0.15]; mk = 'o';
    case 'PV',      col = [0.93 0.60 0.10]; mk = 's';
    case 'Wind',    col = [0.00 0.45 0.74]; mk = '^';
    case 'PV+Wind', col = [0.49 0.18 0.56]; mk = 'd';
    otherwise,      col = [0.5 0.5 0.5];    mk = 'x';
end
end
