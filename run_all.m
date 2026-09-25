%% RUN_ALL - press F5 (Run) to get every result of the project.
%
%  1. Solves the TRANSMISSION system (IEEE 14-bus) and the DISTRIBUTION
%     system (IEEE 33-bus), each with 4 scenarios:
%         Base  |  PV only  |  Wind only  |  PV + Wind
%  2. Prints the result tables in the Command Window.
%  3. Opens all figures (as tabs of one Figures window) and saves them as
%     PNG in results/figures/. Tables are saved as CSV in results/.
%
%  Want the curve of another bus?  Run  explore_buses  afterwards.
%  Want to change plant sizes, buses or settings?  Edit functions/study_config.m
%
%  Author: Hussain Tak (imhussaintak@gmail.com)
%  M.Tech project (2025), EPES, Dept. of Electrical Engineering, NIT Srinagar

clc; close all;
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, 'functions'), fullfile(here, 'data'));
resultsDir = fullfile(here, 'results');
if ~exist(resultsDir, 'dir'), mkdir(resultsDir); end

fprintf('Solving (about 10 s)...\n');
trans = run_study('transmission');
dist  = run_study('distribution');

%% Result tables
show = @(ttl, T) fprintf('\n==== %s ====\n%s', ttl, formattedDisplayText(T, 'SuppressMarkup', true));
show('RENEWABLE PLANTS', [trans.plants; dist.plants]);
show(sprintf('%s  -  renewables at bus %d', upper(trans.cfg.title), trans.cfg.derBus), trans.summary);
show(sprintf('%s  -  renewables at bus %d', upper(dist.cfg.title), dist.cfg.derBus), dist.summary);
fprintf(['LambdaMax = collapse load / base load.  LoadMargin = extra MW the system can take.\n' ...
         'Vmin = lowest load-bus voltage at base load.  Lindex / FVSI close to 1 = near collapse.\n']);

writetable([trans.plants; dist.plants], fullfile(resultsDir, 'renewable_plants.csv'));
writetable(trans.summary, fullfile(resultsDir, 'transmission_summary.csv'));
writetable(dist.summary,  fullfile(resultsDir, 'distribution_summary.csv'));
save(fullfile(resultsDir, 'transmission.mat'), 'trans');
save(fullfile(resultsDir, 'distribution.mat'), 'dist');

%% Figures
make_figures(trans, dist);
fprintf('\nFigures are open (tabs in the Figures window) and saved in %s\n', fullfile(resultsDir, 'figures'));
