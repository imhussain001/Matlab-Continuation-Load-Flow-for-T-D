clc;
clear all;
close all;

% Load data (assumes 'data33' provides bus and branch matrices)
data33;

% Number of buses and branches
bn = length(bus(:,1)); % Number of buses
bnn = length(branch(:,1)); % Number of branches

% Convert data to per-unit (pu) using base voltage 12.66 kV and base MVA = 100
zbase = 12.66 * 12.66 / 100;
branch(:,4) = branch(:,4) ./ zbase; % Convert R to pu
branch(:,5) = branch(:,5) ./ zbase; % Convert X to pu
bus(:,2) = bus(:,2) ./ 100;         % Convert real power (P) to pu
bus(:,3) = bus(:,3) ./ 100;         % Convert reactive power (Q) to pu

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

% Construct BIBC matrix
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
sigma = 0.05; % Step size for continuation
tolerance = 1e-5; % Convergence tolerance
maxIterations = 100; % Max iterations per load flow
maxLambda = 2.0; % Maximum load scaling factor
BusForCPF = 14; % Bus for loadability analysis

% Initialize storage for results
VoltageProfile = [];
LambdaValues = [];

% Continuation Load Flow loop
while lambda <= maxLambda
    % Scale the load with the continuation parameter
    P(:,1) = lambda * complex(bus(:,2), bus(:,3));
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
        % Update voltage using DLF
        % Ensure the correct size for v1
        v1 = DLF(2:bn, 2:bn) * I1; % Voltage drop along branches (size (bn-1) x 1)
        
        % Update bus voltages (excluding slack bus)
        v(2:bn) = v0(2:bn) - v1; % Update bus voltages

        % Update branch currents
        I2 = zeros(bnn, 1);
        for k = 1:bnn
            fromBus = b(k);
            toBus = c(k);
            I2(k) = conj(P(toBus) / v(toBus));
        end

        % Check convergence
        error = max(abs(I1 - I2));
        if error < tolerance
            break; % Convergence achieved
        else
            I1 = I2; % Update for next iteration
        end
    end

    % Store results for plotting
    VoltageProfile = [VoltageProfile; abs(v(BusForCPF))];
    LambdaValues = [LambdaValues; lambda];

    % Break if voltage drops below acceptable limit (e.g., 0.9 pu)
    if abs(v(BusForCPF)) < 0.9
        break;
    end

    % Increment lambda
    lambda = lambda + sigma;
end

% Plot results
figure;
plot(LambdaValues, VoltageProfile, '-o', 'LineWidth', 2);
xlabel('Lambda (Load Scaling Factor)');
ylabel('Voltage Magnitude (pu)');
title(['Loadability Curve for Bus #', num2str(BusForCPF)]);
grid on;
