function Ybus = Calculate_Ybus(BusData, BranchData)
%CALCULATE_YBUS  Bus admittance matrix (pu) from BusData/BranchData.
%
%   Ybus = Calculate_Ybus(BusData, BranchData)
%
%   BusData    : N x 17, columns as documented in data/ieee14_data.m.
%                Uses col 15/16 (shunt G/B, pu).
%   BranchData : L x 21. Uses col 1/2 (from/to), 7/8 (R/X, pu),
%                9 (total line charging B, pu), 15 (off-nominal tap ratio,
%                0 = plain line). The tap is on the FROM side.
%
%   Returns a dense N x N complex matrix (fine for the 14/33-bus systems).
%   Adapted from S. Chatterjee's continuation-power-flow (GPL-2.0).
%
%   Limitation: two parallel branches between the same bus pair are not
%   supported by the off-diagonal loop below.

N = max(BusData(:,1));
Ybus = zeros(N);

% Off-diagonal elements for plain lines (tap == 0)
for i = 1:N
    for j = i+1:N
        pos = (BranchData(:,1)==i & BranchData(:,2)==j) | (BranchData(:,1)==j & BranchData(:,2)==i);
        if any(pos)
            a = BranchData(pos,15);
            if a == 0
                Ybus(i,j) = -1/(BranchData(pos,7) + 1i*BranchData(pos,8));
            end
        end
    end
end
Ybus = Ybus + transpose(Ybus);
Ybus = Ybus + diag(-sum(Ybus,2));

% Line charging (B/2 at each end) and bus shunts
for i = 1:N
    pos = (BranchData(:,1)==i | BranchData(:,2)==i);
    if any(pos)
        Ybus(i,i) = Ybus(i,i) + sum(0.5i .* BranchData(pos,9)) + BusData(i,15) + 1i*BusData(i,16);
    end
end

% Off-nominal tap transformers (tap ratio a on the from side, t = 1/a)
for n = find(BranchData(:,15) > 0).'
    i = BranchData(n,1);
    j = BranchData(n,2);
    t = 1/BranchData(n,15);
    Y = 1./(BranchData(n,7) + 1i*BranchData(n,8));
    Ybus(i,i) = Ybus(i,i) + Y*(abs(t)^2);
    Ybus(i,j) = Ybus(i,j) - conj(t)*Y;
    Ybus(j,i) = Ybus(j,i) - t*Y;
    Ybus(j,j) = Ybus(j,j) + Y;
end
end
