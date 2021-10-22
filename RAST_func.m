function [err, Q_predicted] = RAST_func(x, isotherm, minlnP, M, EoS)
% Real Adsorbed Solution Theory objective function
% x((i-1)*N+j), 1<=j<=N: ln p_j^0 for i-th data point
% x(ndata*N+j): j-th coefficient for EoS
% isotherm(i), 1<=i<=N: function handle for isotherm i, Q(lnP)
% minlnP(i), 1<=i<=N: lnP at which Q is zero
% M{j}(i, 1:2): component j, i-th partial pressure & loading
% EoS: function handle that computes the activity coefficients
%      [\gamma_1, ..., \gamma_N] = EoS([z_1, z_2, ..., z_N])

N = length(M);
ndata = size(M{1}, 1);
lnP = zeros(ndata, N);
Q = zeros(ndata, N);
for i = 1 : N  % components
    lnP(:, i) = log(M{i}(:, 1));
    Q(:, i) = M{i}(:, 2);
end

z = Q ./ sum(Q, 2);
coeff = x(ndata*N+1:end);
lnP0 = reshape(x(1:ndata*N), N, ndata)';
Q_predicted = zeros(ndata, N);
err = 0;
for i = 1 : ndata
    gamma = EoS(coeff, z(i, 1:end-1));
    [y, Q_predicted(i, :)] = RAST_func_per_point([lnP0(i, :), z(i, 1:end-1), gamma], isotherm, minlnP, lnP(i, :), Q(i, :), 2);
    err = err + sum(y.^2);
end

end