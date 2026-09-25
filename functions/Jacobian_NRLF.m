function [WholeJ,J] = Jacobian_NRLF(Ybus,V,d,posSL,posPV)
%JACOBIAN_NRLF  Polar power-flow Jacobian.
%
%   J = [dP/dd   |V|.dP/d|V| ;
%        dQ/dd   |V|.dQ/d|V| ]
%
%   NOTE: the voltage columns are scaled by |V_j| (i.e. derivatives w.r.t.
%   d|V|/|V|). NRLF therefore updates V multiplicatively. Code that wants
%   derivatives w.r.t. |V| itself (e.g. run_cpf) must divide those columns
%   by |V_j|.
%
%   WholeJ is the full 2N x 2N matrix; J has slack-angle rows/cols and
%   slack/PV-magnitude rows/cols removed.

N = length(V);
pos = [posSL; posSL+posPV];

% J11: dPi/d(dj), off-diagonal then diagonal
J11 = -(abs(V)*transpose(abs(V))) .* abs(Ybus) .* ...
    sin(angle(Ybus) + repmat(transpose(d),N,1) - repmat(d,1,N));
J11 = J11.*(1-eye(N));
J11 = J11 - diag(sum(J11,2));

% J21: dQi/d(dj)
J21 = -(abs(V)*transpose(abs(V))) .* abs(Ybus) .* ...
    cos(angle(Ybus) + repmat(transpose(d),N,1) - repmat(d,1,N));
J21 = J21.*(1-eye(N));
J21 = J21 - diag(sum(J21,2));

% J12: |Vj|.dPi/d|Vj|
J12 = J21.*(2*eye(size(J21)) - ones(size(J21)));
J12 = J12 + diag(2*(abs(V).^2).*real(diag(Ybus)));

% J22: |Vj|.dQi/d|Vj|
J22 = J11.*(-2*eye(size(J21)) + ones(size(J21)));
J22 = J22 - diag(2*(abs(V).^2).*imag(diag(Ybus)));

J = [J11 J12; J21 J22];
WholeJ = J;
J(:,pos>0) = [];
J(pos>0,:) = [];
end
