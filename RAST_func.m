function [err, Q_predicted, lnP0, psi] = RAST_func(x, isotherm, minlnP, M, EoS, EoS_deriv, tol)
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

if nargin < 7 || isempty(tol)
    tol = 1e-5;
end

if nargin < 6 || isempty(EoS_deriv)
    % central difference to calculate numerical derivative
    EoS_deriv = @(y)(sum((log(EoS([y(1:end-1), y(end)+tol]))-log(EoS([y(1:end-1), y(end)-tol])))/2/tol.*y(1:end-1)));
end

z = Q ./ sum(Q, 2);
coeff = x(ndata*N+1:end);
EoS_with_coeff = @(y)EoS(coeff, y);
lnP0 = reshape(x(1:ndata*N), N, ndata)';
Q_predicted = zeros(ndata, N);
psi = zeros(ndata, N);
err = 0;
for i = 1 : ndata
    [IAST_err, ~, psi(i, :)] = IAST_func([z(i, 1:end-1), lnP0(i, :)], isotherm, minlnP, lnP(i, :), EoS_with_coeff, [], 2);

    tN = adsorption_potential(isotherm{N}, minlnP(N), lnP0(i, N));
    Q_t = 0;
    for j=1:N
        Q_t = Q_t + z(i, j) / isotherm{j}(lnP0(i, j));
    end
    Q_t = 1 / (Q_t + EoS_deriv(coeff, [z(i, 1:N-1), tN]));
    Q_predicted(i, :) = Q_t * z(i, :);
    
    err = err + sum(IAST_err.^2);
end

err = err + sum(((Q_predicted-Q)./Q).^2, 'all');
end