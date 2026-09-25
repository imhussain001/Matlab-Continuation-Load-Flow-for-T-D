% Define step size and initialize lambda
step_size = 0.01;  
lambda = 0;  
max_steps = 50;  

lambda_values = zeros(max_steps, 1);
voltage_values = zeros(max_steps, 1);
active_power = zeros(max_steps, 1);
reactive_power = zeros(max_steps, 1);

for i = 1:max_steps
    lambda = lambda + step_size; % Increment lambda

    % Call your actual load flow function
    [V, Psch, Qsch] = NRLF();  % Replace this with your solver

    % Extract voltage, active power, and reactive power at Bus 14
    voltage_values(i) = V(14);  
    active_power(i) = Psch(14) - BusData(14, 7);  
    reactive_power(i) = Qsch(14) - BusData(14, 8);

    % Store lambda
    lambda_values(i) = lambda;

    % Stop if voltage collapses
    if V(14) < 0.7
        lambda_values = lambda_values(1:i);
        voltage_values = voltage_values(1:i);
        active_power = active_power(1:i);
        reactive_power = reactive_power(1:i);
        break;
    end
end

% Display table
PowerVoltageTable = table(lambda_values, voltage_values, active_power, reactive_power);
disp(PowerVoltageTable);
