function [err, Q_predicted, x, F] = RAST_func_IAST_solve(x, isotherm, minlnP, M, EoS, mode, x0)
% Real Adsorbed Solution Theory objective function
% x(j): j-th coefficient for EoS
% isotherm(i), 1<=i<=N: function handle for isotherm i, Q(lnP)
% minlnP(i), 1<=i<=N: lnP at which Q is zero
% M{j}(i, 1:2): component j, i-th partial pressure & loading
% EoS: function handle that computes the activity coefficients
%      [\gamma_1, ..., \gamma_N] = EoS([z_1, z_2, ..., z_N])

if nargin < 6 || isempty(mode)
    mode = 2;
end

N = length(M);
ndata = size(M{1}, 1);
P = zeros(ndata, N);
Q = zeros(ndata, N);
for i = 1 : N  % components
    P(:, i) = M{i}(:, 1);
    Q(:, i) = M{i}(:, 2);
end

EoS_with_coeff = @(y)EoS(x, y);

if nargin < 7 || isempty(x0)
    z = zeros(ndata, N-1);
    sumQ = sum(Q, 2);
    for i = 1 : ndata
        z(i, 1:end) = Q(i, 1:end-1)/sumQ(i);
    end
    x0 = IAST_init(log(P), z, EoS_with_coeff(z));
end

[Q_predicted, x, F] = IAST_solve(P, [], isotherm, minlnP, EoS_with_coeff, [], mode, x0);
err = sum(((Q_predicted - Q)./Q).^2, 'all');

end