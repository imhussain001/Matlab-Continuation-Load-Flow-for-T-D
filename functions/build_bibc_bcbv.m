function net = build_bibc_bcbv(fromBus, toBus, z)
%BUILD_BIBC_BCBV  BIBC, BCBV and DLF matrices of a radial feeder (Teng 2003).
%
%   net = build_bibc_bcbv(fromBus, toBus, z)
%
%   fromBus, toBus : L x 1 branch terminals, bus 1 is the substation (slack)
%   z              : L x 1 complex branch impedances (pu)
%
%   BIBC (L x N): I_branch = BIBC * I_injection. BIBC(k,j) = 1 iff bus j
%                 is downstream of branch k.
%   BCBV (N x L): V_1 - V_j = BCBV * I_branch. BCBV(j,k) = z_k iff branch
%                 k lies on the path from the substation to bus j.
%   DLF  (N x N): BCBV * BIBC, so V = V_1 - DLF * I_injection.
%
%   Algorithm (Teng): for branch k = i -> j, column j of BIBC is a copy of
%   column i (every branch feeding i also feeds j) plus a 1 in row k.
%   Branches are processed parent-before-child regardless of input order.
%
%   Legacy bug fixed: the old code accumulated ROWS of BIBC
%   (BIBC(k,:) += BIBC(k-1,:)), which makes each branch carry the current
%   of the upstream buses and ignores laterals. That gave V33 = 0.853 pu
%   at base load instead of the published ~0.913 pu.

fromBus = fromBus(:); toBus = toBus(:); z = z(:);
L = numel(fromBus);
N = max([fromBus; toBus]);
if L ~= N-1
    error('build_bibc_bcbv:NotRadial', ...
        'A radial feeder with %d buses needs %d branches, got %d.', N, N-1, L);
end

BIBC = zeros(L, N);
reached = false(N,1); reached(1) = true;
done = false(L,1);
for pass = 1:L
    progressed = false;
    for k = find(~done).'
        i = fromBus(k); j = toBus(k);
        if reached(i) && ~reached(j)
            BIBC(:,j) = BIBC(:,i);
            BIBC(k,j) = 1;
            reached(j) = true; done(k) = true; progressed = true;
        end
    end
    if all(done), break; end
    if ~progressed
        error('build_bibc_bcbv:NotRadial', ...
            'Branches %s are not connected to the substation as a tree.', mat2str(find(~done).'));
    end
end

net.BIBC = BIBC;
net.BCBV = transpose(BIBC) .* transpose(z);   % BCBV(j,k) = z_k * BIBC(k,j)
net.DLF  = net.BCBV * BIBC;
net.fromBus = fromBus; net.toBus = toBus; net.z = z;
end
