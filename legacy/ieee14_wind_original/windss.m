clc;
clear all;
close all;
% By Hussain Tak (EPES)

%% Step 1: Aerodynamic Model (Steady-State for Grid Integration)

% Wind Turbine Parameters
rho = 1.225;      % Air density (kg/m³)
R = 45;           % Rotor radius (m)
A = pi * R^2;     % Swept area (m²)
Cp_max = 0.4;     % Power coefficient
v_w = 10;         % Wind speed (m/s)

% Mechanical Power from Wind
P_wind_mech = 0.5 * rho * A * Cp_max * v_w^3; % Watts

fprintf('\n------ Aerodynamic Model (Per Turbine) ------\n');
fprintf('Mechanical Power Extracted: %.2f kW\n', P_wind_mech/1000);

%% Step 2: Generator Electrical Output
eta_mech = 0.9; % Mechanical-to-electrical efficiency
P_gen_elec = eta_mech * P_wind_mech; % Active power output (W)

% Assume Reactive Power Generation
Q_gen_elec = 0.3 * P_gen_elec; % VAR (30% assumption)

V_gen = 400; % Generator Output Voltage (V) (initial low voltage)
I_gen = P_gen_elec / V_gen; % Generator Output Current (A)

fprintf('\n------ Generator Model (Per Turbine) ------\n');
fprintf('Electrical Active Power Output: %.2f kW\n', P_gen_elec/1000);
fprintf('Electrical Reactive Power Output: %.2f kVAR\n', Q_gen_elec/1000);
fprintf('Generator Output Voltage: %.2f V\n', V_gen);

%% Step 3: Boost Converter Steady-State Model
Vdc_target = 650; % Target DC Link Voltage
D = 1 - (V_gen / Vdc_target); % Duty Cycle for Boost

% MPPT corrects D if needed (here assume steady-state optimal D)
D = min(max(D, 0.3), 0.6); % Practical limit 30%-60%

Vdc_actual = V_gen / (1 - D); % Boosted DC voltage
I_L = P_gen_elec / V_gen; % Boost inductor current

fprintf('\n------ Boost Converter Model ------\n');
fprintf('Duty Cycle (D): %.2f\n', D);
fprintf('Boosted DC Link Voltage: %.2f V\n', Vdc_actual);

%% Step 4: Inverter and Grid Synchronization
V_grid = 400; % Grid phase voltage (line-to-line)
M = (V_grid * sqrt(2)) / Vdc_actual; % Modulation index

% Ensure Modulation Index limits
M = min(max(M, 0.85), 0.95);

fprintf('\n------ Inverter Model ------\n');
fprintf('Grid Voltage: %.2f V\n', V_grid);
fprintf('Inverter Modulation Index (M): %.2f\n', M);

%% Step 5: System Requirements
P_target = 1.5e6;   % Active Power Required (W) - 1.5 MW
Q_target = 1.8e6;   % Reactive Power Required (VAR) - 1.2 MVAR

% Number of turbines required
N_turbines_P = ceil(P_target / P_gen_elec);
N_turbines_Q = ceil(Q_target / Q_gen_elec);
N_turbines = max(N_turbines_P, N_turbines_Q);

% Total Power Produced
P_total = N_turbines * P_gen_elec;
Q_total = N_turbines * Q_gen_elec;

fprintf('\n------ Wind Farm Total Power ------\n');
fprintf('Number of Turbines: %d\n', N_turbines);
fprintf('Total Active Power Supplied: %.2f MW\n', P_total/1e6);
fprintf('Total Reactive Power Supplied: %.2f MVAR\n', Q_total/1e6);

%% Step 6: Save Data
save('windfarm_with_boost_inverter.mat', 'P_total', 'Q_total', 'N_turbines', 'P_total', 'Q_total', 'Vdc_actual', 'M');

%% Step 7: Plots
figure;
subplot(2,2,1);
bar(P_wind_mech/1000, 'FaceColor', [0.2 0.6 1]);
ylabel('Power (kW)');
title('Mechanical Power per Turbine');
grid on;

subplot(2,2,2);
bar(P_gen_elec/1000, 'FaceColor', [0.8 0.4 0]);
ylabel('Power (kW)');
title('Electrical Power per Turbine');
grid on;

subplot(2,2,3);
bar([P_total/1e6, Q_total/1e6]);
set(gca, 'XTickLabel', {'Total Active Power (MW)', 'Total Reactive Power (MVAR)'});
ylabel('Power (MW/MVAR)');
title('Wind Farm Total Power');
grid on;

subplot(2,2,4);
bar(D, 'FaceColor', [1 0.2 0.2]);
ylabel('Duty Cycle');
title('Boost Converter Duty Cycle');
ylim([0 1]);
grid on;
