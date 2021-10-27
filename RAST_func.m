function [err, Q_predicted, lnP0, psi] = RAST_func(x, isotherm, M, EoS, varargin)
% Real Adsorbed Solution Theory objective function
% x((i-1)*N+j), 1<=j<=N: ln p_j^0 for i-th data point
% x(ndata*N+j): j-th coefficient for EoS
% isotherm(i), 1<=i<=N: function handle for isotherm i, Q(lnP)
% minlnP(i), 1<=i<=N: lnP at which Q is zero
% M{j}(i, 1:2): component j, i-th partial pressure & loading
% EoS: function handle that computes the activity coefficients
%      [\gamma_1, ..., \gamma_N] = EoS([z_1, z_2, ..., z_N])

if nargin < 4 || rem(nargin,2) ~= 0
    error('RAST_func:Number of arguments incorrect');
end

names          = {'minlnP'; 'ads_pot'; 'inv_ads_pot'; 'EoS_deriv'; 'tol'};
default_values = {      [];        [];            [];          []; 1e-5};
opt_args = process_variable_arguments(names, default_values, varargin);
minlnP = opt_args.('minlnP');
ads_pot = opt_args.('ads_pot');
inv_ads_pot = opt_args.('inv_ads_pot');
EoS_deriv = opt_args.('EoS_deriv');
tol = opt_args.('tol');

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
EoS_with_coeff = @(y)EoS(coeff, y);

if isempty(EoS_deriv)
    % central difference to calculate numerical derivative
    EoS_deriv_with_coeff = @(y)(sum((log(EoS_with_coeff([y(1:end-1), y(end)+tol]))-log(EoS_with_coeff([y(1:end-1), y(end)-tol])))/2/tol.*y(1:end-1)));
else
    EoS_deriv_with_coeff = @(y)EoS_deriv(coeff, y);
end

lnP0 = reshape(x(1:ndata*N), N, ndata)';
Q_predicted = zeros(ndata, N);
psi = zeros(ndata, N);
err = 0;
for i = 1 : ndata
    [IAST_err, ~, psi(i, :)] = IAST_func([z(i, 1:end-1), lnP0(i, :)], lnP(i, :), 'mode', 1, 'isotherm', isotherm, 'minlnP', minlnP, 'tol', tol, 'EoS', EoS_with_coeff, 'ads_pot', ads_pot, 'inv_ads_pot', inv_ads_pot);

    Q_t = 0;
    for j=1:N
        Q_t = Q_t + z(i, j) / isotherm{j}(lnP0(i, j));
    end
    Q_t = 1 / (Q_t + EoS_deriv_with_coeff([z(i, 1:N-1), psi(i, N)]));
    Q_predicted(i, :) = Q_t * z(i, :);
    
    err = err + sum(IAST_err.^2);
end

err = err + sum(((Q_predicted-Q)./Q).^2, 'all');
end