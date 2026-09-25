clc;
clear;
close all;

%% System Data (33-bus radial distribution system)
% Note: Enter the line data (from bus, to bus, impedance) and load data (P, Q) here.

% Example placeholder data (replace with actual system data)
line_data = [
    1, 2, 0.0922 + 0.0470j;
    2, 3, 0.4930 + 0.2511j;
    3, 4, 0.3660 + 0.1864j;
    4, 5, 0.3811 + 0.1941j;
    % Add more lines as needed
];
load_data = [
    2, 100, 60; % Bus, P (kW), Q (kVAR)
    3, 90, 40;
    4, 120, 80;
    5, 60, 30;
    % Add more loads as needed
];

n_bus = max(max(line_data(:, 1:2))); % Number of buses
n_lines = size(line_data, 1);

% Initialize system parameters
base_power = 100; % Base power in kVA
base_voltage = 12.66; % Base voltage in kV
slack_bus = 1; % Slack bus index

% Initialize load multiplier
lambda_max = 2.5; % Maximum load multiplier for continuation
delta_lambda = 0.05; % Load increment step

%% Backward/Forward Sweep Method
V_base = base_voltage * ones(n_bus, 1); % Initialize voltage
V = V_base; % Initial voltage at all buses
tolerance = 1e-4; % Convergence tolerance

% Initialize results
lambda_values = 0:delta_lambda:lambda_max;
voltages = zeros(n_bus, length(lambda_values));

% Main loop for continuation load flow
for idx = 1:length(lambda_values)
    lambda = lambda_values(idx);
    
    % Update load demands
    P_load = zeros(n_bus, 1);
    Q_load = zeros(n_bus, 1);
    for i = 1:size(load_data, 1)
        bus = load_data(i, 1);
        P_load(bus) = lambda * load_data(i, 2) / base_power; % Per-unit P
        Q_load(bus) = lambda * load_data(i, 3) / base_power; % Per-unit Q
    end
    
    % Iterative Backward/Forward Sweep
    converged = false;
    while ~converged
        % Backward Sweep (Calculate branch currents)
        I_line = zeros(n_lines, 1);
        I_bus = zeros(n_bus, 1);
        for l = n_lines:-1:1
            from = line_data(l, 1);
            to = line_data(l, 2);
            z = line_data(l, 3);
            
            % Update currents
            I_line(l) = conj(P_load(to) + 1j * Q_load(to)) / V(to) + I_bus(to);
            I_bus(from) = I_bus(from) + I_line(l);
        end
        
        % Forward Sweep (Update voltages)
        V_new = V_base;
        for l = 1:n_lines
            from = line_data(l, 1);
            to = line_data(l, 2);
            z = line_data(l, 3);
            
            % Update voltages
            V_new(to) = V_new(from) - z * I_line(l);
        end
        
        % Check convergence
        if max(abs(V_new - V)) < tolerance
            converged = true;
        end
        V = V_new;
    end
    
    % Store results
    voltages(:, idx) = abs(V);
end

%% Plot Loadability Curve
figure;
for bus = 2:n_bus % Exclude slack bus
    plot(lambda_values, voltages(bus, :), 'LineWidth', 1.5);
    hold on;
end
xlabel('Load Multiplication Factor (\lambda)');
ylabel('Voltage Magnitude (pu)');
title('Loadability Curve for 33-bus System');
legend(arrayfun(@(x) sprintf('Bus %d', x), 2:n_bus, 'UniformOutput', false));
grid on;
