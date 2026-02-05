% Clear workspace and close figures
clear all;
close all;

%% Step 1: Solar PV System Modeling

% Constants
K = 1.38065e-23; % Boltzmann Constant (J/K)
q = 1.602e-19;   % Electron charge (C)
Iscn = 8.21;     % Short-circuit current (A)
Vocn = 32.9;     % Open-circuit voltage (V)
Ns = 72;         % Number of series cells per panel
T = 25 + 273.15; % Operating temperature in Kelvin
Gn = 1000;       % Irradiance at STC (W/m^2)
G = 1000;        % Actual irradiance (W/m^2)
Rs = 0.15;       % Series resistance (Ohms)
Rp = 500;        % Parallel resistance (Ohms)

% Thermal voltage
Vtn = Ns * (K * T) / q;

% Reverse saturation current
I0 = Iscn / (exp(Vocn / (2 * Vtn)) - 1);

% Photocurrent
Ipv = (G / Gn) * Iscn;

% Generate I-V curve using Newton-Raphson method
V = linspace(0, Vocn, 100);
I = zeros(size(V));

for k = 1:length(V)
    Ik = Ipv; % Initial guess for current
    for j = 1:10 % Newton-Raphson iterations
        f = Ik - Ipv + I0 * (exp((V(k) + Ik * Rs) / (Vtn * 2)) - 1) + (V(k) + Ik * Rs) / Rp;
        df = 1 + I0 * (Rs / (Vtn * 2)) * exp((V(k) + Ik * Rs) / (Vtn * 2)) + Rs / Rp;
        Ik = Ik - f / df;
    end
    I(k) = Ik;
end

% Power calculation
P = V .* I;
[~, idx] = max(P);
Vmp_single = V(idx); % Voltage at MPP for single panel
Ipv_mpp = I(idx);    % Current at MPP for single panel
Ppv_mpp = P(idx);    % Power at MPP for single panel

%% Step 2: Scaling to 1 MW

% Desired total active power
P_target = 1.5e6; % 1 MW

% Calculate total number of panels required
N_total = ceil(P_target / Ppv_mpp);

% Desired system voltage (e.g., for inverter input)
Vpv_target = 440; % Adjust as needed

% Calculate number of panels in series
N_series = ceil(Vpv_target / Vmp_single);

% Calculate number of parallel strings
N_parallel = ceil(N_total / N_series);

% Update total number of panels
N_total = N_series * N_parallel;

% Total array voltage and current
Vpv = N_series * Vmp_single;
Ipv_total = N_parallel * Ipv_mpp;

% Total active power
P_pv = Vpv * Ipv_total;

%% Step 3: Reactive Power Calculation

% Desired power factor
pf = 0.8071; % For equal active and reactive power (cos(45°))

% Apparent power
S_pv = P_pv / pf;

% Reactive power
Q_pv = sqrt(S_pv^2 - P_pv^2);

%% Step 4: Boost Converter Modeling

% Target DC voltage
Vdc_target = 650; % Target output voltage
D = 1 - (Vpv / Vdc_target);
D = min(max(D, 0.3), 0.6); % Limit duty cycle
Vdc_actual = Vpv / (1 - D); % Boosted voltage

%% Step 5: Inverter Conversion and Grid Synchronization

V_grid = 400; % Grid voltage
M = (V_grid * sqrt(2)) / Vdc_actual;
M = min(max(M, 0.85), 0.95); % Ensure modulation index is within range

%% Step 6: Power Exchange with Grid

P_steady_state = P_pv;  % PV supplies active power
Q_steady_state = Q_pv;  % PV supplies reactive power

% Save data
save('pv_data.mat', 'P_pv', 'Q_pv');

%% Step 7: Plot Results in 2x2 Subplots

figure;

% Subplot 1 (Top-Left): I-V Characteristics
subplot(2,2,1);
plot(V, I, 'r', 'LineWidth', 2); hold on;
plot(Vmp_single, Ipv_mpp, 'bo', 'MarkerFaceColor', 'b');
xlabel('Voltage (V)');
ylabel('Current (A)');
title('I-V Characteristics of Solar PV');
grid on;

% Subplot 2 (Top-Right): P-V Characteristics
subplot(2,2,2);
plot(V, P, 'k', 'LineWidth', 2); hold on;
plot(Vmp_single, Ppv_mpp, 'go', 'MarkerFaceColor', 'g');
xlabel('Voltage (V)');
ylabel('Power (W)');
title('P-V Characteristics of Solar PV');
grid on;

% Subplot 3 (Bottom-Left): Power Bar Chart (Active & Reactive)
subplot(2,2,3);
bar([P_steady_state, Q_steady_state] / 1e6); % Convert to MW/MVAR
set(gca, 'XTickLabel', {'Active Power (MW)', 'Reactive Power (MVAR)'});
ylabel('Power (MW/MVAR)');
title('Total Power Supplied by PV System');
grid on;

% Subplot 4 (Bottom-Right): Duty Cycle of Boost Converter
subplot(2,2,4);
bar(D, 'FaceColor', [1 0.2 0.2]);
ylabel('Duty Cycle');
title('Boost Converter Duty Cycle');
ylim([0 1]);
grid on;


%% Step 8: Display Results

fprintf('\n------ PV Array Configuration ------\n');
fprintf('Total Panels: %d\n', N_total);
fprintf('Series Panels: %d\n', N_series);
fprintf('Parallel Strings: %d\n', N_parallel);
fprintf('Total Array Voltage: %.2f V\n', Vpv);
fprintf('Total Array Current: %.2f A\n', Ipv_total);
fprintf('Total Active Power Supplied by PV: %.2f MW\n', P_pv/1e6);
fprintf('Total Reactive Power Supplied by PV: %.2f MVAR\n', Q_pv/1e6);
fprintf('Boost Converter Output Voltage: %.2f V\n', Vdc_actual);
fprintf('Inverter Modulation Index: %.2f\n', M);
