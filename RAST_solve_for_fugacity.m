function [err, gamma, lnP0, psi] = RAST_solve_for_fugacity(x, M, EoS, varargin)
% Real Adsorbed Solution Theory objective function
% x(j): j-th coefficient for EoS
% isotherm(i), 1<=i<=N: function handle for isotherm i, Q(lnP)
% minlnP(i), 1<=i<=N: lnP at which Q is zero
% M{j}(i, 1:2): component j, i-th partial pressure & loading
% EoS: function handle that computes the activity coefficients
%      [\gamma_1, ..., \gamma_N] = EoS([z_1, z_2, ..., z_N])
% mode: 6 minimizes the squared errors in ln(P) and S_1j
%       7 minimizes the squared errors in ln(P_i)

if nargin < 3 || rem(nargin,2) ~= 1
    error('RAST_solve_for_fugacity:Number of arguments incorrect');
end

options_0 = optimset('FinDiffType', 'central', 'FunValCheck', 'on', 'MaxFunEvals', 800, 'MaxIter', 100, 'TolFun', 1e-6, 'TolX', 1e-6, 'Display', 'off');
names          = {'isotherm'; 'minlnP'; 'ads_pot'; 'inv_ads_pot'; 'EoS_deriv'; 'tol'; 'options'; 'mode'};
default_values = {        [];       [];        [];            [];          [];  1e-5; options_0; 6};
opt_args = process_variable_arguments(names, default_values, varargin);
isotherm = opt_args.('isotherm');
minlnP = opt_args.('minlnP');
ads_pot = opt_args.('ads_pot');
inv_ads_pot = opt_args.('inv_ads_pot');
EoS_deriv = opt_args.('EoS_deriv');
tol = opt_args.('tol');
options = opt_args.('options');
mode = opt_args.('mode');

N = length(M);
ndata = size(M{1}, 1);
lnP = zeros(ndata, N);
Q = zeros(ndata, N);

if isempty(isotherm)
    isotherm = cell(1, N);
end
if isempty(ads_pot)
    ads_pot = cell(1, N);
end
if isempty(inv_ads_pot)
    inv_ads_pot = cell(1, N);
end

for i = 1 : N  % components
    lnP(:, i) = log(M{i}(:, 1));
    Q(:, i) = M{i}(:, 2);
    if isempty(inv_ads_pot{i})
        if isempty(ads_pot{i})
            if isempty(isotherm{i}) || isempty(minlnP)
                error('RAST_solve_for_fugacity:isotherm and minlnP must be provided without inv_ads_pot');
            end
            ads_pot{i} = @(y)adsorption_potential(y, isotherm{i}, minlnP(i), tol);
        end
        inv_ads_pot{i} = @(y)inv_adsorption_potential(y, ads_pot{i});
    end
end

Q_t = sum(Q, 2);
z = Q ./ Q_t;
options = optimset(options_0, options);
EoS_with_coeff = @(y)EoS(x, y);

if isempty(EoS_deriv)
    % central difference to calculate numerical derivative
    EoS_deriv_with_coeff = @(y)(sum((log(EoS_with_coeff([y(1:end-1), y(end)+tol]))-log(EoS_with_coeff([y(1:end-1), y(end)-tol])))/2/tol.*y(1:end-1)));
else
    EoS_deriv_with_coeff = @(y)EoS_deriv(x, y);
end

function err = excess_loading(y, Q_t, z, isotherm, inv_ads_pot, EoS_deriv_with_coeff)
    % x: adsorption potential
    ex_inv_Q_t = 1 / Q_t;
    for k = 1 : length(isotherm)
        ex_inv_Q_t = ex_inv_Q_t - z(k) / isotherm{k}(inv_ads_pot{k}(y));
    end
    % ex_inv_Q_t
    % EoS_deriv_with_coeff([z(1:end-1), y])
    err = EoS_deriv_with_coeff([z(1:end-1), y]) - ex_inv_Q_t;
end

psi = zeros(ndata, 1);
gamma = zeros(ndata, N);
lnP0 = zeros(ndata, N);
lnP_p = zeros(ndata, N);
for i = 1 : ndata
    func = @(y)excess_loading(y, Q_t(i), z(i, :), isotherm, inv_ads_pot, EoS_deriv_with_coeff);
    x0 = 0;
    [psi(i), resnorm, residual, exitflag, output, lambda, jacobian] = lsqnonlin(func, x0, 0, [], options);
    gamma(i, :) = EoS_with_coeff([z(i, 1:end-1), psi(i)]);
    for j = 1 : N
        lnP0(i, j) = inv_ads_pot{j}(psi(i));
    end
    lnP_p(i, :) = lnP0(i, :) + log(gamma(i, :) .* z(i, :));
end

if mode == 6
    P_p = exp(lnP_p);
    P_p_t = sum(P_p, 2);
    y_p = P_p ./ P_p_t;
    selectivity_p = (z(:, 1) ./ y_p(:, 1)) ./ (z(:, 2:end) ./ y_p(:, 2:end));
    P = exp(lnP);
    P_t = sum(P, 2);
    y = P ./ P_t;
    selectivity = (z(:, 1) ./ y(:, 1)) ./ (z(:, 2:end) ./ y(:, 2:end));
    err_P = sum((log(P_p_t) - log(P_t)).^2);
    err_S = sum((selectivity_p - selectivity).^2, 'all');
    err = err_P + err_S;
elseif mode == 7
    err = sum((lnP_p - lnP).^2, 'all');
else
    error('RAST_solve_for_fugacity:UnknownMode','Mode parameter outside known choices.')
end
end