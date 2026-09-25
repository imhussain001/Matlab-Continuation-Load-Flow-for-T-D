function pb = cpf_power_balance(cpf, BusData, Sder, Sbase_MVA)
%CPF_POWER_BALANCE  System power balance at every CPF point (MW / Mvar).
%
%   pb = cpf_power_balance(cpf, BusData, Sder, Sbase_MVA)
%
%   cpf    : result of run_cpf
%   Sder   : N x 1 constant complex DER injection (pu, generation positive)
%
%   Fields (M x 1, one value per CPF point, in MW / Mvar):
%     lambda, loadP, loadQ        total system load  lambda*(PL, QL)
%     derP, derQ                  DER injection (constant)
%     slackP, slackQ              slack generator / substation supply
%                                 (negative = reverse flow to the upstream grid)
%     genP, genQ                  all generators incl. slack (Q from cpf.Qg)
%     lossP, lossQ                network losses (sum of all injections)
%     derShare                    derP / loadP in % (Inf at no load)
%     upper                       logical, points from lambda = 0 up to the nose

M = numel(cpf.lambda);
lam = cpf.lambda(:);
Vc = cpf.V .* exp(1i*cpf.delta);            % N x M
Sinj = Vc .* conj(cpf.Ybus*Vc);             % net injections, N x M (pu)
slack = find(cpf.posSL, 1);

PL = BusData(:,7); QL = BusData(:,8);
pb.lambda = lam;
pb.loadP = lam*sum(PL)*Sbase_MVA;
pb.loadQ = lam*sum(QL)*Sbase_MVA;
pb.derP = repmat(sum(real(Sder))*Sbase_MVA, M, 1);
pb.derQ = repmat(sum(imag(Sder))*Sbase_MVA, M, 1);
pb.slackP = (real(Sinj(slack,:)).' + lam*PL(slack) - real(Sder(slack)))*Sbase_MVA;
pb.slackQ = cpf.Qg(slack,:).'*Sbase_MVA;
gen = cpf.posPV | cpf.posSL;
pb.genP = (sum(real(Sinj(gen,:)),1).' + lam*sum(PL(gen)) - sum(real(Sder(gen))))*Sbase_MVA;
pb.genQ = sum(cpf.Qg(gen,:),1).'*Sbase_MVA;
pb.lossP = sum(real(Sinj),1).'*Sbase_MVA;
pb.lossQ = sum(imag(Sinj),1).'*Sbase_MVA;
pb.derShare = 100*pb.derP./pb.loadP;
pb.upper = (1:M).' <= cpf.idxMax;
end
