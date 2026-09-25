function [L, loadBuses] = l_index(Ybus, Vc, busType)
%L_INDEX  Kessel-Glavitsch voltage stability L-index of each load bus.
%
%   [L, loadBuses] = l_index(Ybus, Vc, busType)
%
%   Ybus    : N x N bus admittance matrix (pu) - use the SAME Ybus as the
%             load flow (taps and shunts included).
%   Vc      : N x 1 complex bus voltages from a solved load flow (pu).
%   busType : N x 1, 3 = slack, 2 = PV, 0 = PQ.
%
%       F   = -Y_LL^{-1} Y_LG
%       L_j = | 1 - sum_i F_ji V_i / V_j |,   j in load buses
%
%   L -> 0 at no load and L -> 1 approaching voltage collapse.
%
%   Fixes w.r.t. legacy ind.m: uses the solved operating point instead of
%   the textbook base-case voltages, keeps the complex sum inside |1 - .|
%   (legacy took |1 - |sum||), and reuses Calculate_Ybus so transformer
%   taps and the bus-9 shunt are included.

isGen = busType == 3 | busType == 2;
loadBuses = find(~isGen);
genBuses  = find(isGen);
F = -Ybus(loadBuses,loadBuses) \ Ybus(loadBuses,genBuses);
L = abs(1 - (F*Vc(genBuses)) ./ Vc(loadBuses));
end
