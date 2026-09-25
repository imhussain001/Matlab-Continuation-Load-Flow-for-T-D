function res = bfs_pv_curve(net, Sload, Sfixed, opts)
%BFS_PV_CURVE  Upper branch of the P-V curve by repeated BFS load flows.
%
%   res = bfs_pv_curve(net, Sload, Sfixed, Step=0.1, MinStep=1e-4)
%
%   Solves  Snet(lambda) = lambda*Sload - Sfixed  for lambda = 0, Step, ...
%   warm-starting each solve from the previous one. When a solve fails the
%   step is halved (down to MinStep), so the last point approaches the
%   loadability limit from below.
%
%   This is a *repeated* load flow, not a continuation method: BFS has no
%   way to pass the nose, so only the stable (upper) branch is obtained and
%   lambdaMax is a lower bound. Use run_cpf for the full nose curve.
%
%   Sload  : N x 1 complex base-case load (pu), load positive
%   Sfixed : N x 1 complex constant DER injection (pu), generation positive
%
%   Output: lambda (M x 1), V (N x M complex), lambdaMax (last converged).

arguments
    net struct
    Sload double
    Sfixed double
    opts.Step (1,1) double {mustBePositive} = 0.1
    opts.MinStep (1,1) double {mustBePositive} = 1e-4
    opts.MaxPoints (1,1) double {mustBeInteger, mustBePositive} = 10000
    opts.Tolerance (1,1) double {mustBePositive} = 1e-10
    opts.MaxIterations (1,1) double {mustBeInteger, mustBePositive} = 500
end

N = size(net.DLF,1);
lam = zeros(opts.MaxPoints,1);
Vh  = complex(nan(N, opts.MaxPoints));

[Vc, ok] = bfs_load_flow(net, -Sfixed, Tolerance=opts.Tolerance, MaxIterations=opts.MaxIterations);
if ~ok
    error('bfs_pv_curve:NoStartingPoint', 'BFS did not converge at lambda = 0.');
end
n = 1; lam(1) = 0; Vh(:,1) = Vc;

step = opts.Step;
while step >= opts.MinStep && n < opts.MaxPoints
    lamTry = lam(n) + step;
    [Vtry, ok] = bfs_load_flow(net, lamTry*Sload - Sfixed, VInit=Vh(:,n), ...
        Tolerance=opts.Tolerance, MaxIterations=opts.MaxIterations);
    if ok
        n = n + 1; lam(n) = lamTry; Vh(:,n) = Vtry;
    else
        step = step/2;
    end
end

res.lambda = lam(1:n);
res.V = Vh(:,1:n);
res.lambdaMax = lam(n);
end
