clear all;
close all;

% Data for 33-Bus Radial Distribution System
branch = [1  1  2  0.0922 0.0470;
          2  2  3  0.4930 0.2511;
          3  3  4  0.3660 0.1864;
          4  4  5  0.3811 0.1941;
          5  5  6  0.8190 0.7070;
          6  6  7  0.1872 0.6188;
          7  7  8  0.7114 0.2351;
          8  8  9  1.0300 0.7400;
          9  9  10 1.0040 0.7400;
          10 10 11 0.1996 0.0650;
          11 11 12 0.3744 0.1238;
          12 12 13 1.4680 1.1550;
          13 13 14 0.5416 0.7129;
          14 14 15 0.5910 0.5260;
          15 15 16 0.7463 0.5450;
          16 16 17 1.2890 1.7210;
          17 17 18 0.7320 0.5740;
          18 2  19 0.1640 0.1565;
          19 19 20 1.5042 1.3554;
          20 20 21 0.4095 0.4784;
          21 21 22 0.7089 0.9373;
          22 3  23 0.4512 0.3083;
          23 23 24 0.8980 0.7091;
          24 24 25 0.8960 0.7011;
          25 6  26 0.2030 0.1034;
          26 26 27 0.2842 0.1447;
          27 27 28 1.0590 0.9337;
          28 28 29 0.8042 0.7006;
          29 29 30 0.5075 0.2585;
          30 30 31 0.9744 0.9630;
          31 31 32 0.3105 0.3619;
          32 32 33 0.3410 0.5302];
      
bus = [1  0.000  0.000;
       2  0.100  0.060;
       3  0.090  0.040;
       4  0.120  0.080;
       5  0.060  0.030;
       6  0.060  0.020;
       7  0.200  0.100;
       8  0.200  0.100;
       9  0.060  0.020;
       10 0.060  0.020;
       11 0.045  0.030;
       12 0.060  0.035;
       13 0.060  0.035;
       14 0.120  0.080;
       15 0.060  0.010;
       16 0.060  0.020;
       17 0.060  0.020;
       18 0.090  0.040;
       19 0.090  0.040;
       20 0.090  0.040;
       21 0.090  0.040;
       22 0.090  0.040;
       23 0.090  0.050;
       24 0.420  0.200;
       25 0.420  0.200;
       26 0.060  0.025;
       27 0.060  0.025;
       28 0.060  0.020;
       29 0.120  0.070;
       30 0.200  0.600;
       31 0.150  0.070;
       32 0.210  0.100;
       33 0.060  0.040];

% Parameters
reactivePowerValues = [0.2, 0.2, 0.6, 0.8]; % Reactive power values to test at Bus 33 (in MW)
colors = ['r', 'g', 'b', 'm']; % Colors for different curves
BusForCPF = 33; % Bus for loadability analysis
sigma = 0.1; % Step size for continuation
tolerance = 1e-5; % Convergence tolerance
maxLambda = 3.2; % Maximum load scaling factor
maxIterations = 100; % Max iterations per load flow

% Convert base voltage and power to per-unit
zbase = 12.66 * 12.66 / 100;

% Loop through each reactive power value
figure;
hold on;
for idx = 1:length(reactivePowerValues)
    % Update reactive power at Bus 33
    bus(33, 3) = reactivePowerValues(idx); % Set Q value
    bus(:, 2:3) = bus(:, 2:3) / 100; % Convert to per-unit
    
    % Convert branch impedance to per-unit
    branch(:,4:5) = branch(:,4:5) / zbase;

    % Initialize variables
    b = branch(:,2); % From bus
    c = branch(:,3); % To bus
    v = complex(ones(length(bus),1), 0); % Flat start
    BIBC = zeros(length(branch), length(bus));
    BCBV = zeros(length(bus), length(branch));
    z = complex(branch(:,4), branch(:,5));
    
    % Construct BIBC matrix
    for k = 1:length(branch)
        i = b(k);
        j = c(k);
        BIBC(k, j) = 1;
        if k > 1
            BIBC(k, :) = BIBC(k, :) + BIBC(k-1, :);
        end
    end
    
    % Construct BCBV matrix
    for k = 1:length(branch)
        i = b(k);
        j = c(k);
        BCBV(j, k) = z(k);
        if i ~= 0
            BCBV(j, :) = BCBV(i, :) + BCBV(j, :);
        end
    end
    
    % Calculate DLF matrix
    DLF = BCBV * BIBC;
    
    % Initialize continuation parameters
    lambda = 0;
    VoltageAtBus = [];
    LambdaValues = [];
    
    % Continuation Load Flow Loop
    while lambda <= maxLambda
        % Scale reactive power only
        P = lambda * complex(bus(:,2), bus(:,3));
        
        % Iterative Backward/Forward Sweep
        I1 = zeros(length(branch), 1);
        for iter = 1:maxIterations
            % Update branch currents
            for k = 1:length(branch)
                fromBus = b(k);
                toBus = c(k);
                I1(k) = conj(P(toBus) / v(toBus));
            end
            
            % Voltage drop calculation
            v1 = DLF(2:end, 2:end) * I1;
            v(2:end) = v(1) - v1; % Update voltage (excluding slack bus)
            
            % Check for convergence
            if max(abs(I1 - conj(P(2:end) ./ v(2:end)))) < tolerance
                break;
            end
        end
        
        % Store results
        VoltageAtBus = [VoltageAtBus; abs(v(BusForCPF))];
        LambdaValues = [LambdaValues; lambda];
        
        lambda = lambda + sigma; % Increment lambda
    end
    
    % Plot the loadability curve for current reactive power value
    plot(LambdaValues, VoltageAtBus, [colors(idx), 'o-'], 'LineWidth', 1.5, 'DisplayName', ...
         ['Q = ', num2str(reactivePowerValues(idx)), ' MW']);
end

% Customize plot
xlabel('Lambda (Load Scaling Factor)');
ylabel('Voltage Magnitude at Bus 33 (pu)');
title('Impact of Reactive Power on Voltage Stability at Bus 33');
grid on;
legend('show');
hold off;
