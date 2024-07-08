function [err, Q_predicted, x_IAST, err_IAST, lnP0_IAST, psi_IAST] = RAST_func_IAST_solve(x, isotherm, M, EoS, varargin)
% Real Adsorbed Solution Theory objective function
% x(j): j-th coefficient for EoS
% isotherm(i), 1<=i<=N: function handle for isotherm i, Q(lnP)
% minlnP(i), 1<=i<=N: lnP at which Q is zero
% M{j}(i, 1:2): component j, i-th partial pressure & loading
% EoS: function handle that computes the activity coefficients
%      [\gamma_1, ..., \gamma_N] = EoS([z_1, z_2, ..., z_N-1, \Psi])
% EoS_deriv: function handle that computes the excess inverse loading
%            d(G^ex/RT) / d\Psi
%   (1/Q_t)^excess = EoS_deriv([coeff_1, ..., coeff_M], [z_1, ..., z_N-1, \Psi])

if nargin < 4 || rem(nargin,2) ~= 0
    error('RAST_func:Number of arguments incorrect');
end

options_0 = optimset('Display', 'off');
names          = {'minlnP'; 'maxlnP'; 'ads_pot'; 'inv_ads_pot'; 'EoS_deriv'; 'mode'; 'x0'; 'x_lb'; 'x_ub'; 'tol'; 'options'};
default_values = {      [];       [];        [];            [];          [];      1;   [];     [];     [];  1e-5; options_0};
opt_args = process_variable_arguments(names, default_values, varargin);
minlnP = opt_args.('minlnP');
maxlnP = opt_args.('maxlnP');
ads_pot = opt_args.('ads_pot');
inv_ads_pot = opt_args.('inv_ads_pot');
EoS_deriv = opt_args.('EoS_deriv');
mode = opt_args.('mode');
x0 = opt_args.('x0');
x_lb = opt_args.('x_lb');
x_ub = opt_args.('x_ub');
tol = opt_args.('tol');
options = opt_args.('options');

N = length(M);
ndata = size(M{1}, 1);
P = zeros(ndata, N);
Q = zeros(ndata, N);
for i = 1 : N  % components
    P(:, i) = M{i}(:, 1);
    Q(:, i) = M{i}(:, 2);
end

options = optimset(options_0, options);
EoS_with_coeff = @(y)EoS(x, y);

if isempty(EoS_deriv)
    % central difference to calculate numerical derivative
    EoS_deriv_with_coeff = @(y)(sum((log(EoS_with_coeff([y(1:end-1), y(end)+tol]))-log(EoS_with_coeff([y(1:end-1), y(end)-tol])))/2/tol.*[y(1:end-1),1-sum(y(1:end-1))]));
else
    EoS_deriv_with_coeff = @(y)EoS_deriv(x, y);
end

if isempty(x0)
    z = zeros(ndata, N-1);
    sumQ = sum(Q, 2);
    for i = 1 : ndata
        z(i, 1:end) = Q(i, 1:end-1)/sumQ(i);
    end
    x0 = IAST_init(log(P), z, EoS_with_coeff([z,1e3*ones(ndata,1)]), mode);
end

[Q_predicted, x_IAST, err_IAST, lnP0_IAST, psi_IAST, exitflags_IAST] = IAST_solve(P, [], 'isotherm', isotherm, 'minlnP', minlnP, 'maxlnP', maxlnP, 'EoS', EoS_with_coeff, 'ads_pot', ads_pot, 'inv_ads_pot', inv_ads_pot, 'options', options, 'tol', tol, 'mode', mode, 'x0', x0, 'x_lb', x_lb, 'x_ub', x_ub, 'EoS_deriv', EoS_deriv_with_coeff);
if ~strcmp(optimget(options, 'Display'), 'off')
    disp('exitflags_IAST: ')
    disp(exitflags_IAST)
end
diff_sq = (relative_error_safe(Q_predicted, Q)).^2;
diff_sq_valid = diff_sq(exitflags_IAST > 0, :);
err = sum(diff_sq_valid, 'all') / numel(diff_sq_valid);

end