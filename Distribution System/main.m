clear all;
close all;

 
% Impedance matrix (complex R + jX)
z(:,1) = complex(branch(:,4), branch(:,5));

% Define bus connection parameters
b = branch(:,2); % From bus
c = branch(:,3); % To bus

% Initialize voltage vector (complex) with magnitude 1 and angle 0
v(1:bn,1) = complex(1,0);
v0 = v; % Flat start

% Program for BIBC (Bus Injection to Branch Current) and BCBV matrices
BIBC = zeros(bnn, bn); % bnn = number of branches, bn = number of buses
BCBV = zeros(bn, bnn); % bn = number of buses, bnn = number of branches


%% Load power values
load('pv_data.mat', 'P_pv', 'Q_pv');  % Load saved PV power values
load('windfarm_with_boost_inverter.mat', 'P_total', 'Q_total');  % Load saved Wind power values 

%% Convert  power to per-unit
base_MVA = 100e6;  
P_pv_pu = (P_pv + P_total) / base_MVA;  
Q_pv_pu = (Q_pv + Q_total) / base_MVA;

% Inject PV power at Bus 33
bus(33,2) = bus(33,2) - P_pv_pu; % Reduce demand by PV active power
bus(33,3) = bus(33,3) - Q_pv_pu; % Reduce demand by PV reactive power

fprintf('PV and Wind System added at Bus 33: P = %.2f MW, Q = %.2f MVAR\n', (P_pv + P_total)/1e6, (Q_pv + Q_total)/1e6);
%% Construct BIBC matrix
for k = 1:bnn
    i = b(k); % From bus
    j = c(k); % To bus
    BIBC(k, j) = 1; % Update current injection relation
    if k > 1
        BIBC(k, :) = BIBC(k, :) + BIBC(k - 1, :); % Update for k > 1
    end
end

% Construct BCBV matrix
for k = 1:bnn
    i = b(k); % From bus
    j = c(k); % To bus
    BCBV(j, k) = z(k); % Update branch impedance
    if i ~= 0
        BCBV(j, :) = BCBV(i, :) + BCBV(j, :);
    end
end

% Calculate DLF (Distribution Load Flow) matrix
DLF = BCBV * BIBC; % (bn x bnn) * (bnn x bn) -> (bn x bn)

% Initialize continuation parameters
lambda = 0; % Starting continuation parameter
sigma = 0.1; % Step size for continuation
tolerance = 1e-5; % Convergence tolerance
maxIterations = 100; % Max iterations per load flow
maxLambda = 3.5; % Maximum load scaling factor
BusForCPF = 33; % Bus for loadability analysis

% Initialize storage for results
VoltageAtBus = []; % Voltage magnitudes at the selected bus
LambdaValues = []; % Lambda values

%% Continuation Load Flow loop
while lambda <= maxLambda
    % Scale the load with the continuation parameter
    P(:,1) = lambda * complex(bus(:,2), bus(:,3)); % P and Q scaling
    P1 = P(2:bn); % Exclude slack bus

    % Initialize current injection for branches
    I1 = zeros(bnn, 1);
    for k = 1:bnn
        fromBus = b(k);
        toBus = c(k);
        I1(k) = conj(P(toBus) / v(toBus)); % Current in each branch
    end

    % Iterative Backward/Forward Sweep
    for iter = 1:maxIterations
        % Update voltage using DLF considering both real and reactive power
        v1 = DLF(2:bn, 2:bn) * I1; % Voltage drop along branches (size (bn-1) x 1)
        
        % Update bus voltages (excluding slack bus)
        v(2:bn) = v0(2:bn) - v1; % Update bus voltages

        % Update branch currents based on new voltage values
        I2 = zeros(bnn, 1);
        for k = 1:bnn
            fromBus = b(k);
            toBus = c(k);
            I2(k) = conj(P(toBus) / v(toBus));
        end

        % Convergence check
        if max(abs(I1 - I2)) < tolerance
            break; % Converged
        end

        I1 = I2; % Update for next iteration
    end

    % Store results at each step
    VoltageAtBus = [VoltageAtBus; abs(v(BusForCPF))]; % Voltage magnitude at chosen bus
    LambdaValues = [LambdaValues; lambda]; % Scaling factor for load

    lambda = lambda + sigma; % Increase lambda for next iteration
end

%% Plot Loadability Curve
figure;
plot(LambdaValues, VoltageAtBus, 'go-', 'LineWidth', 1.5, 'MarkerSize', 6); % Voltage vs Lambda
hold on;
plot(LambdaValues, 1 - VoltageAtBus, 'r^-', 'LineWidth', 1.5, 'MarkerSize', 6); % Complementary curve
xlabel('Lambda (Load Scaling Factor)');
ylabel('Voltage Magnitude (pu)');
title(['Loadability Curve for Bus #', num2str(BusForCPF)]);
grid on;
legend('Voltage Profile', 'Complementary Curve');

% Store voltage profile before and after PV injection
VoltageWithPV = abs(v); % Voltage after PV injection

