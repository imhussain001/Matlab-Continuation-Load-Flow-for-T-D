function [Vc, converged, iter, Ibr] = bfs_load_flow(net, Snet, opts)
%BFS_LOAD_FLOW  Backward/forward sweep load flow using the DLF matrix.
%
%   [Vc, converged, iter, Ibr] = bfs_load_flow(net, Snet, VInit=..., ...)
%
%   net  : struct from build_bibc_bcbv
%   Snet : N x 1 complex net LOAD at each bus (pu), load positive,
%          DER generation negative. Entry 1 (substation) is ignored.
%
%   Iterates  I = conj(S./V),  V = V_slack - DLF*I  until max|dV| < Tol.
%   Constant-power loads. Returns converged = false instead of silently
%   storing a non-converged voltage (legacy behaviour past the nose).
%
%   Name-value: VInit (N x 1, default flat 1.0), Vslack (1.0),
%               Tolerance (1e-10 pu), MaxIterations (500)

arguments
    net struct
    Snet double
    opts.VInit double = []
    opts.Vslack (1,1) double = 1
    opts.Tolerance (1,1) double {mustBePositive} = 1e-10
    opts.MaxIterations (1,1) double {mustBeInteger, mustBePositive} = 500
end

N = size(net.DLF,1);
Snet = Snet(:);
Vc = opts.VInit;
if isempty(Vc), Vc = complex(ones(N,1)); end
Vc = Vc(:); Vc(1) = opts.Vslack;

converged = false;
for iter = 1:opts.MaxIterations
    I = conj(Snet./Vc);
    I(1) = 0;
    Vnew = opts.Vslack - net.DLF*I;
    if ~all(isfinite(Vnew)), break; end
    dV = max(abs(Vnew - Vc));
    Vc = Vnew;
    if dV < opts.Tolerance
        converged = true;
        break
    end
end
I = conj(Snet./Vc); I(1) = 0;
Ibr = net.BIBC * I;
end
