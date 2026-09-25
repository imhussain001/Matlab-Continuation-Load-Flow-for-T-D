function FVSI = fvsi_index(fromBus, toBus, z, Vc, Ibr)
%FVSI_INDEX  Fast Voltage Stability Index (Musirin & Abdul Rahman, 2002).
%
%   FVSI = fvsi_index(fromBus, toBus, z, Vc, Ibr)
%
%   fromBus, toBus : L x 1 branch terminals
%   z              : L x 1 complex branch impedance (pu)
%   Vc             : N x 1 complex bus voltages from a SOLVED load flow (pu)
%   Ibr            : L x 1 complex branch currents, from -> to (pu)
%
%       FVSI_k = 4 |Z_k|^2 Q_r,k / ( |V_s,k|^2 X_k )
%
%   Q_r is the reactive power arriving at the receiving end of branch k
%   (the whole downstream demand + losses, not just the local bus load).
%   FVSI -> 1 indicates the line is at its voltage-stability limit;
%   FVSI < 0 means reactive power is flowing back towards the source.
%
%   Fixes w.r.t. legacy indccds.m: legacy used the receiving-bus LOAD Q
%   and a flat 1.0 pu sending voltage, which ignores all downstream flow,
%   DER injections and loading level.

Sr = Vc(toBus) .* conj(Ibr);          % complex power into the receiving bus
Qr = imag(Sr);
Vs = abs(Vc(fromBus));
FVSI = 4*abs(z).^2 .* Qr ./ (Vs.^2 .* imag(z));
end
