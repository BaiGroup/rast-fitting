function [c, ceq] = RAST_const(x, isotherm, minlnP, M, EoS)
% Real Adsorbed Solution Theory constraint function
% x((i-1)*N+j), 1<=j<=N: ln p_j^0 for i-th data point
% x(ndata*N+(i-1)*(N-1)+j), 1<=j<=N-1: z_j for component j in the absorbed phase
% x(ndata*(2N-1)+j): j-th coefficient for EoS
% isotherm(i), 1<=i<=N: function handle for isotherm i, Q(lnP)
% minlnP(i), 1<=i<=N: lnP at which Q is zero
% M{j}(i, 1:2): component j, i-th partial pressure & loading
% EoS: function handle that computes the activity coefficients
%      [\gamma_1, ..., \gamma_N] = EoS([z_1, z_2, ..., z_N])

N = length(M);
ndata = size(M{1}, 1);
lnP = zeros(ndata, N);
for i = 1 : N  % components
    lnP(:, i) = log(M{i}(:, 1));
end

coeff = x(ndata*(2*N-1)+1:end);
lnP0 = reshape(x(1:ndata*N), N, ndata)';
z = reshape(x(ndata*N+1:ndata*(2*N-1)), N-1, ndata)';
err = zeros(ndata, 2*N-1);
for i = 1 : ndata
    gamma = EoS(coeff, z(i, :));
    err(i, :) = IAST_func([lnP0(i, :), z(i, :)], lnP(i, :), 'isotherm', isotherm, 'minlnP', minlnP, 'EoS', gamma, 'mode', 2);
end

c = [];
ceq = reshape(err, ndata*(2*N-1), 1);

end