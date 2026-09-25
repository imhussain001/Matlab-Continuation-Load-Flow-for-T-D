function make_figures(trans, dist, opts)
%MAKE_FIGURES  All result figures, shown on screen and saved as PNG.
%
%   make_figures(trans, dist)                  trans/dist from run_study
%   make_figures(trans, dist, Show=false)      save only (no windows)
%
%   Every figure compares the four scenarios Base / PV / Wind / PV+Wind
%   with the same colours. Figures open as tabs of one window (docked);
%   use Dock=false for separate windows. PNGs go to results/figures/.
%
%    1  Transmission - loadability (P-V) curve
%    2  Transmission - voltage profile at base load
%    3  Transmission - losses and slack generation vs load
%    4  Transmission - L-index (voltage stability index)
%    5  Transmission - generator reactive limits
%    6  Distribution - loadability (P-V) curves
%    7  Distribution - voltage profile at base load
%    8  Distribution - losses and substation power vs load
%    9  Distribution - FVSI (voltage stability index)
%   10  Renewable power injected (P and Q of each plant)
%   11  Loadability summary (lambda_max and load margin)
%   12  PV and wind plant models

arguments
    trans struct
    dist struct
    opts.Show (1,1) logical = true
    opts.Dock (1,1) logical = true
    opts.Save (1,1) logical = true
end

figDir = fullfile(project_root(), 'results', 'figures');
if opts.Save && ~exist(figDir, 'dir'), mkdir(figDir); end
ctx = struct('show', opts.Show, 'dock', opts.Dock, 'save', opts.Save, 'dir', figDir, 'n', 0);

ctx = fig_loadability(ctx, trans);
ctx = fig_voltage_profile(ctx, trans);
ctx = fig_losses(ctx, trans);
ctx = fig_index(ctx, trans);
ctx = fig_qlimits(ctx, trans);
ctx = fig_loadability(ctx, dist);
ctx = fig_voltage_profile(ctx, dist);
ctx = fig_losses(ctx, dist);
ctx = fig_index(ctx, dist);
ctx = fig_injection(ctx, trans, dist);
ctx = fig_summary(ctx, trans, dist);
fig_models(ctx, trans, dist);
end

% ========================================================================
function ctx = fig_loadability(ctx, out)
cfg = out.cfg; sc = out.scenarios;
buses = cfg.plotBuses;
[ctx, f] = new_figure(ctx, sprintf('%s - loadability curve', cfg.label), [100 100 560*numel(buses)+80 470]);
t = tiledlayout(f, 1, numel(buses), 'TileSpacing', 'compact');
for b = buses
    ax = nexttile(t); hold(ax, 'on');
    for k = 1:numel(sc)
        c = sc(k).cpf; [col, mk] = style(sc(k).name);
        plot(ax, c.lambda, c.V(b,:), '-', 'Color', col, 'LineWidth', 1.8, ...
            'DisplayName', sprintf('%s  (\\lambda_{max} = %.3f)', sc(k).name, c.lambdaMax));
        plot(ax, c.lambdaMax, c.V(b,c.idxMax), mk, 'Color', col, 'MarkerFaceColor', col, ...
            'MarkerSize', 7, 'HandleVisibility', 'off');
    end
    xline(ax, 1, ':', 'base load', 'HandleVisibility', 'off', 'LabelVerticalAlignment', 'bottom');
    what = ternary(b == cfg.derBus, 'renewable bus', 'weakest bus');
    decorate(ax, 'Load multiplier \lambda  (1 = base load)', sprintf('Voltage at bus %d (p.u.)', b), ...
        sprintf('Bus %d (%s)', b, what));
    ylim(ax, [0 inf]);
end
title(t, sprintf('%s: loadability (P-V) curves   (marker = collapse point)', cfg.title));
ctx = finish(ctx, f, sprintf('%s_loadability', lower(cfg.label)));
end

function ctx = fig_voltage_profile(ctx, out)
cfg = out.cfg; sc = out.scenarios;
[ctx, f] = new_figure(ctx, sprintf('%s - voltage profile', cfg.label));
ax = axes(f); hold(ax, 'on');
N = numel(sc(1).op.V);
for k = 1:numel(sc)
    [col, mk] = style(sc(k).name);
    plot(ax, 1:N, sc(k).op.V, ['-' mk], 'Color', col, 'MarkerFaceColor', col, 'MarkerSize', 4, ...
        'LineWidth', 1.4, 'DisplayName', sc(k).name);
end
yline(ax, 0.95, '--', '0.95 p.u.', 'HandleVisibility', 'off');
yline(ax, 1.05, '--', '1.05 p.u.', 'HandleVisibility', 'off');
decorate(ax, 'Bus number', 'Voltage (p.u.)', sprintf('%s: bus voltages at base load (\\lambda = 1)', cfg.title));
xlim(ax, [1 N]);
ctx = finish(ctx, f, sprintf('%s_voltage_profile', lower(cfg.label)));
end

function ctx = fig_losses(ctx, out)
cfg = out.cfg; sc = out.scenarios; sys = out.system;
source = ternary(strcmp(cfg.system, 'ieee33'), 'Substation', 'Slack generator');
[ctx, f] = new_figure(ctx, sprintf('%s - losses and grid supply', cfg.label), [100 100 1100 440]);
t = tiledlayout(f, 1, 2, 'TileSpacing', 'compact');
ax1 = nexttile(t); hold(ax1, 'on');
ax2 = nexttile(t); hold(ax2, 'on');
for k = 1:numel(sc)
    Sder = complex(zeros(size(sys.BusData,1),1));
    Sder(cfg.derBus) = complex(sc(k).der.P_MW, sc(k).der.Q_Mvar)/sys.Sbase_MVA;
    pb = cpf_power_balance(sc(k).cpf, sys.BusData, Sder, sys.Sbase_MVA);
    u = pb.upper;
    col = style(sc(k).name);
    plot(ax1, pb.lambda(u), pb.lossP(u), 'Color', col, 'LineWidth', 1.6, 'DisplayName', sc(k).name);
    plot(ax2, pb.lambda(u), pb.slackP(u), 'Color', col, 'LineWidth', 1.6, 'DisplayName', sc(k).name);
end
xline(ax1, 1, ':', 'HandleVisibility', 'off'); xline(ax2, 1, ':', 'HandleVisibility', 'off');
yline(ax2, 0, 'k-', 'HandleVisibility', 'off');
decorate(ax1, 'Load multiplier \lambda', 'Active power losses (MW)', 'Network losses');
decorate(ax2, 'Load multiplier \lambda', sprintf('%s output (MW)', source), ...
    sprintf('%s active power (below 0 = power sent back)', source));
title(t, sprintf('%s: losses and grid supply from no load to collapse', cfg.title));
ctx = finish(ctx, f, sprintf('%s_losses', lower(cfg.label)));
end

function ctx = fig_index(ctx, out)
cfg = out.cfg; sc = out.scenarios;
idx = sc(1).index;
vals = cell2mat(arrayfun(@(s) s.index.value(:), sc, 'UniformOutput', false).');
[ctx, f] = new_figure(ctx, sprintf('%s - %s', cfg.label, idx.name), [100 100 1000 460]);
ax = axes(f);
b = bar(ax, categorical(idx.labels), vals, 'grouped');
for k = 1:numel(sc), b(k).FaceColor = style(sc(k).name); end
if strcmp(idx.name, 'FVSI')
    xl = 'Branch (labelled by its receiving-end bus)';
    note = 'FVSI \rightarrow 1 = line at its stability limit;  < 0 = reactive power flowing back';
else
    xl = 'Load bus';
    note = 'L \rightarrow 1 = voltage collapse';
end
decorate(ax, xl, idx.name, sprintf('%s: %s at base load   (%s)', cfg.title, idx.name, note));
legend(ax, cfg.scenarios, 'Location', 'southoutside', 'Orientation', 'horizontal');
ctx = finish(ctx, f, sprintf('%s_%s', lower(cfg.label), lower(strrep(idx.name, '-', ''))));
end

function ctx = fig_qlimits(ctx, out)
sc = out.scenarios; cfg = out.cfg;
if isempty(sc(1).cpfNoQLim), return; end
c = sc(1).cpf; c0 = sc(1).cpfNoQLim; bus = cfg.plotBuses(1);
[ctx, f] = new_figure(ctx, sprintf('%s - generator Q limits', cfg.label), [100 100 1100 440]);
t = tiledlayout(f, 1, 2, 'TileSpacing', 'compact');
ax = nexttile(t); hold(ax, 'on');
plot(ax, c0.lambda, c0.V(bus,:), '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.4, ...
    'DisplayName', sprintf('Unlimited generator Q  (\\lambda_{max} = %.3f)', c0.lambdaMax));
plot(ax, c.lambda, c.V(bus,:), 'k-', 'LineWidth', 1.8, ...
    'DisplayName', sprintf('With Q limits  (\\lambda_{max} = %.3f)', c.lambdaMax));
hit = c.events(startsWith(c.events.Action, "PV->PQ"), :);
plot(ax, c.lambda(hit.Point), c.V(bus, hit.Point), 'rv', 'MarkerFaceColor', 'r', ...
    'DisplayName', 'A generator reaches Q_{max}');
decorate(ax, 'Load multiplier \lambda', sprintf('Voltage at bus %d (p.u.)', bus), 'Base case with and without Q limits');
ylim(ax, [0 inf]);
ax = nexttile(t); hold(ax, 'on');
gens = find(c.posPV); gc = lines(numel(gens));
up = (1:numel(c.lambda)).' <= c.idxMax;
for g = 1:numel(gens)
    b = gens(g);
    plot(ax, c.lambda(up), out.system.Sbase_MVA*c.Qg(b,up), '-', 'Color', gc(g,:), 'LineWidth', 1.5, ...
        'DisplayName', sprintf('Generator at bus %d', b));
    yline(ax, out.system.Sbase_MVA*c.qLimits.QgMax(b), ':', 'Color', gc(g,:), 'LineWidth', 1.2, 'HandleVisibility', 'off');
end
decorate(ax, 'Load multiplier \lambda', 'Reactive power output (Mvar)', 'Generator Q up to collapse (dotted = Q_{max})');
title(t, sprintf('%s: generator reactive power limits (Base case)', cfg.title));
ctx = finish(ctx, f, sprintf('%s_generator_q_limits', lower(cfg.label)));
end

function ctx = fig_injection(ctx, trans, dist)
[ctx, f] = new_figure(ctx, 'Renewable power injected', [100 100 1100 720]);
t = tiledlayout(f, 2, 2, 'TileSpacing', 'compact');
systems = {trans, dist};
for s = 1:2
    out = systems{s}; sc = out.scenarios;
    P = [arrayfun(@(x) pw(x.der.pv, 'P_W'), sc), arrayfun(@(x) pw(x.der.wind, 'P_W'), sc)]/1e6;
    Q = [arrayfun(@(x) pw(x.der.pv, 'Q_W'), sc), arrayfun(@(x) pw(x.der.wind, 'Q_W'), sc)]/1e6;
    baseLoad = sum(out.system.BusData(:,7))*out.system.Sbase_MVA;
    data = {P, Q}; unit = {'MW', 'Mvar'}; what = {'Active', 'Reactive'};
    for m = 1:2
        ax = nexttile(t);
        cats = categorical(out.cfg.scenarios, out.cfg.scenarios);
        b = bar(ax, cats, data{m}, 'stacked');
        b(1).FaceColor = style('PV'); b(2).FaceColor = style('Wind');
        tot = sum(data{m}, 2);
        if m == 1
            lbl = compose('%.2f MW\n(%.1f%% of load)', tot, 100*tot/baseLoad);
        else
            lbl = compose('%.2f Mvar', tot);
        end
        text(ax, 1:numel(tot), tot, lbl, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 9);
        ylim(ax, [0 max(tot)*1.3]);
        grid(ax, 'on'); box(ax, 'on');
        ylabel(ax, sprintf('%s power (%s)', what{m}, unit{m}));
        title(ax, sprintf('%s: %s power at bus %d', out.cfg.label, lower(what{m}), out.cfg.derBus), 'FontWeight', 'normal');
        legend(ax, {'PV plant', 'Wind farm'}, 'Location', 'northwest');
    end
end
title(t, 'Renewable power injected by each plant (constant, independent of load)');
ctx = finish(ctx, f, 'renewable_injection');
end

function ctx = fig_summary(ctx, trans, dist)
[ctx, f] = new_figure(ctx, 'Loadability summary', [100 100 1100 450]);
t = tiledlayout(f, 1, 2, 'TileSpacing', 'compact');
systems = {trans, dist};
for s = 1:2
    out = systems{s}; sc = out.scenarios;
    lm = arrayfun(@(x) x.cpf.lambdaMax, sc);
    ax = nexttile(t);
    b = bar(ax, categorical(out.cfg.scenarios, out.cfg.scenarios), lm, 'FaceColor', 'flat');
    for k = 1:numel(sc), b.CData(k,:) = style(sc(k).name); end
    gain = 100*(lm/lm(1) - 1);
    text(ax, 1:numel(lm), lm, compose('%.3f\n(%+.1f%%)', lm, gain), 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', 'FontSize', 9);
    ylim(ax, [0 max(lm)*1.25]);
    grid(ax, 'on'); box(ax, 'on');
    ylabel(ax, '\lambda_{max}  (collapse load / base load)');
    title(ax, out.cfg.title, 'FontWeight', 'normal');
end
title(t, 'Maximum loadability \lambda_{max}   (label: change vs Base)');
ctx = finish(ctx, f, 'loadability_summary');
end

function ctx = fig_models(ctx, trans, dist)
pv = first_plant(trans, 'pv');
[ctx, f] = new_figure(ctx, 'PV and wind plant models', [100 100 1100 640]);
t = tiledlayout(f, 2, 2, 'TileSpacing', 'compact');
ax = nexttile(t); plot(ax, pv.curve.V, pv.curve.I, 'LineWidth', 1.6, 'Color', style('PV')); hold(ax, 'on');
plot(ax, pv.Vmp_V, pv.Imp_A, 'ko', 'MarkerFaceColor', 'k');
grid(ax, 'on'); xlabel(ax, 'Panel voltage (V)'); ylabel(ax, 'Panel current (A)');
title(ax, 'PV panel I-V curve (1000 W/m^2, 25 °C)', 'FontWeight', 'normal');
ax = nexttile(t); plot(ax, pv.curve.V, pv.curve.P, 'LineWidth', 1.6, 'Color', style('PV')); hold(ax, 'on');
plot(ax, pv.Vmp_V, pv.Pmp_W, 'ko', 'MarkerFaceColor', 'k');
grid(ax, 'on'); xlabel(ax, 'Panel voltage (V)'); ylabel(ax, 'Panel power (W)');
title(ax, sprintf('PV panel P-V curve (maximum power point %.0f W)', pv.Pmp_W), 'FontWeight', 'normal');
ax = nexttile(t, [1 2]); axis(ax, 'off');
lines_ = strings(0,1);
for out = {trans, dist}
    o = out{1}; p = first_plant(o, 'pv'); w = first_plant(o, 'wind');
    lines_(end+1) = sprintf('%s (plants at bus %d):', o.cfg.title, o.cfg.derBus); %#ok<AGROW>
    lines_(end+1) = sprintf('    PV plant : %d panels (%d in series x %d strings) = %.2f MW + j%.2f Mvar, pf %.3f', ...
        p.N_panels, p.N_series, p.N_parallel, p.P_W/1e6, p.Q_W/1e6, p.PowerFactor); %#ok<AGROW>
    lines_(end+1) = sprintf('    Wind farm: %d turbines x %.2f MW = %.2f MW + j%.2f Mvar, pf %.3f', ...
        w.N_turbines, w.PerTurbineP_W/1e6, w.P_W/1e6, w.Q_W/1e6, w.P_W/hypot(w.P_W, w.Q_W)); %#ok<AGROW>
    lines_(end+1) = ""; %#ok<AGROW>
end
text(ax, 0, 0.5, lines_, 'FontSize', 11, 'VerticalAlignment', 'middle', 'FontName', 'Consolas');
title(t, 'Renewable plant models');
finish(ctx, f, 'plant_models');
end

% ========================================================================
function [col, mk] = style(scenario)
[col, mk] = scenario_style(scenario);
end

function [ctx, f] = new_figure(ctx, name, pos)
if nargin < 3, pos = [100 100 820 500]; end
ctx.n = ctx.n + 1;
f = figure('Color', 'w', 'Position', pos, 'Name', sprintf('%d. %s', ctx.n, name), ...
    'NumberTitle', 'off', 'Visible', ternary(ctx.show, 'on', 'off'));
if ctx.show && ctx.dock
    f.WindowStyle = 'docked';
end
end

function decorate(ax, xl, yl, ttl)
grid(ax, 'on'); box(ax, 'on');
xlabel(ax, xl); ylabel(ax, yl); title(ax, ttl, 'FontWeight', 'normal');
set(ax, 'FontSize', 10);
legend(ax, 'Location', 'best');
end

function ctx = finish(ctx, f, slug)
if ctx.save
    save_png(f, fullfile(ctx.dir, sprintf('%02d_%s.png', ctx.n, slug)));
end
if ~ctx.show
    close(f);
end
end

function p = first_plant(out, kind)
sc = out.scenarios;
k = find(arrayfun(@(s) ~isempty(s.der.(kind)), sc), 1);
p = sc(k).der.(kind);
end

function v = pw(plant, field)
if isempty(plant), v = 0; else, v = plant.(field); end
end

function v = ternary(c, a, b)
if c, v = a; else, v = b; end
end
