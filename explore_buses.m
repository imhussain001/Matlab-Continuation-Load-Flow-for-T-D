function explore_buses(systemName, buses, opts)
%EXPLORE_BUSES  Loadability (P-V) curves for any buses you choose.
%
%   explore_buses                              asks everything in the Command
%                                              Window; press Enter to stop
%   explore_buses('distribution')              asks only for the buses
%   explore_buses('transmission', [4 9 14])    no questions
%   explore_buses('distribution', 18, XAxis="MW", Scenarios=["Base","PV+Wind"], Save=true)
%
%   Bus input accepts:  14  |  4 9 14  |  10-14, 18  |  all  |  pq
%
%   Each selected bus gets its own panel with the P-V curves of all
%   scenarios, the nose point, the base-load operating point and the
%   0.95 p.u. limit. A table with |V| at base load and at the nose is
%   printed to the command window. No load flows are re-solved: this
%   reads results/transmission.mat or results/distribution.mat written by
%   run_all (it offers to solve the system first if the file is missing).
%
%   Saved figures go to results/figures/buses/.
%   InputFcn replaces MATLAB's input() (e.g. to script the prompts in a test).
%
%   Author: Hussain Tak (imhussaintak@gmail.com)
%   M.Tech project (2025), EPES, Dept. of Electrical Engineering, NIT Srinagar

arguments
    systemName (1,1) string = ""
    buses double = []
    opts.XAxis (1,1) string {mustBeMember(opts.XAxis, ["", "lambda", "MW"])} = ""
    opts.Scenarios string = strings(0)
    opts.Save logical {mustBeScalarOrEmpty} = logical.empty
    opts.InputFcn function_handle {mustBeScalarOrEmpty} = function_handle.empty   % for scripted tests
end

if isempty(opts.InputFcn)
    opts.InputFcn = @(prompt) input(prompt, 's');
end
root = fileparts(mfilename('fullpath'));
addpath(fullfile(root, 'functions'), fullfile(root, 'data'));
interactive = isempty(buses);

if systemName == ""
    choice = ask(opts.InputFcn, 'Which system?  [1] Transmission (IEEE 14-bus)   [2] Distribution (IEEE 33-bus)  [1]: ', "1");
    systemName = ternary(choice == "2", "distribution", "transmission");
end
out = load_study(char(systemName), opts.InputFcn);
if isempty(out), return; end
B = out.system.BusData;

% Ask only for options that were not passed in
if opts.XAxis == ""
    if interactive
        ax = ask(opts.InputFcn, 'X-axis:  [1] loading factor lambda   [2] total system load (MW)  [1]: ', "1");
        opts.XAxis = ternary(ax == "2", "MW", "lambda");
    else
        opts.XAxis = "lambda";
    end
end
if isempty(opts.Save)
    if interactive
        s = ask(opts.InputFcn, 'Save figures as PNG?  (y/N): ', "n");
        opts.Save = startsWith(lower(s), "y");
    else
        opts.Save = false;
    end
end

keep = true(1, numel(out.scenarios));
if ~isempty(opts.Scenarios)
    keep = ismember(lower(string({out.scenarios.name})), lower(opts.Scenarios));
    if ~any(keep)
        error('explore_buses:NoScenario', 'None of %s are scenarios (%s).', ...
            strjoin(opts.Scenarios, ', '), strjoin(string({out.scenarios.name}), ', '));
    end
end
sc = out.scenarios(keep);

while true
    if interactive
        gens = find(B(:,4)==2).';
        fprintf('\n%s: buses 1-%d (slack: %d, generators: %s, renewables at bus %d)\n', out.cfg.title, ...
            size(B,1), find(B(:,4)==3), ternary(isempty(gens), 'none', num2str(gens)), out.cfg.derBus);
        str = ask(opts.InputFcn, 'Buses to plot (e.g. 14 | 4 9 14 | 10-14 | all | pq), Enter to quit: ', "");
        if str == "", break; end
        try
            buses = parse_bus_list(str, B(:,4));
        catch err
            fprintf(2, '%s\n', err.message);
            continue
        end
    end
    plot_buses(out, sc, buses, opts);
    print_table(out, sc, buses);
    if ~interactive, break; end
end
end

% ------------------------------------------------------------------------
function out = load_study(systemName, inputFcn)
cfg = study_config(systemName);
label = lower(cfg.label);                 % 'transmission' / 'distribution'
f = fullfile(project_root(), 'results', [label '.mat']);
if isfile(f)
    S = load(f);
    out = S.(ternary(strcmp(label, 'transmission'), 'trans', 'dist'));
    return
end
s = ask(inputFcn, sprintf('No saved results for the %s system. Solve it now (about 5 s)?  (Y/n): ', label), "y");
if startsWith(lower(s), "n")
    out = [];
    return
end
out = run_study(label);
end

function plot_buses(out, sc, buses, opts)
B = out.system.BusData;
baseLoad_MW = sum(B(:,7))*out.system.Sbase_MVA;
typeName = containers.Map({3, 2, 0}, {'slack', 'generator', 'load'});

perFig = 6;
for first = 1:perFig:numel(buses)
    chunk = buses(first:min(first+perFig-1, numel(buses)));
    n = numel(chunk);
    cols = min(n, 3); rows = ceil(n/cols);
    f = figure('Color', 'w', 'Position', [80 80 380*cols+60 330*rows+60], ...
        'Name', sprintf('%s - buses %s', out.cfg.label, mat2str(chunk)), 'NumberTitle', 'off');
    t = tiledlayout(f, rows, cols, 'TileSpacing', 'compact', 'Padding', 'compact');
    for b = chunk
        ax = nexttile(t); hold(ax, 'on');
        for k = 1:numel(sc)
            c = sc(k).cpf; [col, mk] = scenario_style(sc(k).name);
            x = xdata(c.lambda, opts.XAxis, baseLoad_MW);
            plot(ax, x, c.V(b,:), '-', 'Color', col, 'LineWidth', 1.4, 'DisplayName', sc(k).name);
            plot(ax, x(c.idxMax), c.V(b,c.idxMax), mk, 'Color', col, ...
                'MarkerFaceColor', col, 'HandleVisibility', 'off');
            plot(ax, xdata(1, opts.XAxis, baseLoad_MW), sc(k).op.V(b), 'x', 'Color', col, ...
                'MarkerSize', 8, 'LineWidth', 1.5, 'HandleVisibility', 'off');
        end
        yline(ax, 0.95, 'k--', 'HandleVisibility', 'off');
        xline(ax, xdata(1, opts.XAxis, baseLoad_MW), 'k:', 'HandleVisibility', 'off');
        grid(ax, 'on'); box(ax, 'on'); ylim(ax, [0 inf]);
        xlabel(ax, ternary(opts.XAxis == "MW", 'Total system load (MW)', 'Load multiplier \lambda  (1 = base load)'));
        ylabel(ax, sprintf('Voltage at bus %d (p.u.)', b));
        title(ax, sprintf('Bus %d (%s bus)', b, typeName(B(b,4))), 'FontWeight', 'normal');
    end
    lg = legend(nexttile(t, 1), 'Location', 'southwest');
    lg.FontSize = 8;
    title(t, sprintf('%s: loadability curves  (marker = collapse point, x = base load, dashed = 0.95 p.u.)', out.cfg.title));
    if opts.Save
        d = fullfile(project_root(), 'results', 'figures', 'buses');
        if ~exist(d, 'dir'), mkdir(d); end
        name = sprintf('%s_buses_%s_%s.png', lower(out.cfg.label), strjoin(string(chunk), '-'), opts.XAxis);
        save_png(f, fullfile(d, name), 300);
        fprintf('Saved %s\n', fullfile(d, name));
    end
end
end

function print_table(out, sc, buses)
rows = numel(buses)*numel(sc);
Bus = zeros(rows,1); Scenario = strings(rows,1);
V_base = zeros(rows,1); V_nose = zeros(rows,1); LambdaMax = zeros(rows,1); dV_pct = zeros(rows,1);
r = 0;
for b = buses
    for k = 1:numel(sc)
        r = r + 1;
        Bus(r) = b; Scenario(r) = sc(k).name;
        V_base(r) = sc(k).op.V(b);
        V_nose(r) = sc(k).cpf.V(b, sc(k).cpf.idxMax);
        LambdaMax(r) = sc(k).cpf.lambdaMax;
        dV_pct(r) = 100*(V_base(r) - sc(1).op.V(b))/sc(1).op.V(b);
    end
end
T = table(Bus, Scenario, V_base, V_nose, LambdaMax, dV_pct);
T.Properties.VariableNames{'dV_pct'} = sprintf('dV_vs_%s_pct', matlab.lang.makeValidName(sc(1).name));
fprintf('\n');
disp(T);
end

function x = xdata(lambda, mode, baseLoad_MW)
if mode == "MW", x = lambda*baseLoad_MW; else, x = lambda; end
end

function s = ask(inputFcn, prompt, default)
s = strtrim(string(inputFcn(prompt)));
if s == "", s = default; end
end

function v = ternary(c, a, b)
if c, v = a; else, v = b; end
end
