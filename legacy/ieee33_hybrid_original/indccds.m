% Assume branch and bus data already loaded and converted to per-unit

Vs = ones(bn,1);  % Approximate voltage magnitude (flat start or from load flow)

FVSI = zeros(bnn,1);
for i = 1:bnn
    from = branch(i,2);  % sending bus
    to   = branch(i,3);  % receiving bus
    R = branch(i,4);
    X = branch(i,5);
    Z = sqrt(R^2 + X^2);
    Q_to = bus(to,3);       % reactive load at receiving bus
    V_s = Vs(from);         % sending voltage (assumed 1.0 pu for simplicity)

    FVSI(i) = (4 * Z^2 * Q_to) / (V_s^2 * X);
end

% Plot FVSI
figure;
bar(1:bnn, FVSI, 'FaceColor', [0.4 0.6 0.8]);
grid on;
xlabel('Branch Index');
ylabel('FVSI');
title('Fast Voltage Stability Index (FVSI) for 33-Bus Distribution System');
ylim([0, max(0.7, max(FVSI)+0.1)]);
yticks(0:0.2:ceil(max(FVSI)));
set(gca, 'FontSize', 12);

% Highlight branches with FVSI > 1
hold on;
critical_idx = find(FVSI > 0.7);
bar(critical_idx, FVSI(critical_idx), 'r');
legend('Stable Branches', 'Critical (FVSI > 0.7)', 'Location', 'northwest');

% Display critical branches
disp('?? Critical Branches (FVSI > 0.7):');
disp(table(critical_idx, FVSI(critical_idx), ...
     'VariableNames', {'BranchNumber', 'FVSI'}));
% === Create FVSI Summary Table ===
branch_nums = (1:bnn)';                  % Branch indices
from_buses = branch(:,2);                % Sending buses
to_buses = branch(:,3);                  % Receiving buses
is_critical = FVSI > 0.7;                  % Logical flag

% Build table
FVSI_Table = table(branch_nums, from_buses, to_buses, FVSI, is_critical, ...
    'VariableNames', {'Branch', 'From', 'To', 'FVSI', 'Critical'});

% Display first few rows
disp('--- FVSI Summary Table ---');
disp(FVSI_Table);

