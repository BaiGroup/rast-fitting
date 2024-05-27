% IAST_solve Ideal Adsorbed Solution Theory (IAST) modelling
%
% Syntax
%   Q = IAST_solve(M, S)
%   Q = IAST_solve(M, S, 'tol', 1e-5, 'options', options)
%
% Description
%   For a R-component system (R>1),
%   M is a P x R matrix and M(i,:) is a mixture state point, characterized
%     by R values of partial pressures of the corresponding component.
%   S is a cell array of length R, and S{i} stores the single-component
%     isotherm data for the i-th component. S{i} is a T x 2 matrix, where
%     the first column contains pressures and the second column loadings.
%     T can be different for different components.
%   tol is passed to IAST_func and defaults to 1e-5. It controls the
%     precision of integrals via the trapezoidal method.
%   options are passed to lsqnonlin and can be contructed using optimset.
%   isotherm(i), 1<=i<=N: function handle for isotherm i, Q(lnP)
%   minlnP(i), 1<=i<=N: lnP at which Q is zero
%   EoS: function handle that computes the activity coefficients
%     [\gamma_1, ..., \gamma_N] = EoS([z_1, z_2, ..., z_N-1, \Psi])
%   EoS_deriv: function handle that computes the excess inverse loading
%              d(G^ex/RT) / d\Psi
%     (1/Q_t)^excess = EoS_deriv([coeff_1, ..., coeff_M], [z_1, ..., z_N-1, \Psi])
%   mode: 1 or -1 uses fsolve() and P*y1/P_i^0 - x_i*gamma_i = 0
%     2 or -2 uses lsqnonlin() and ln(P*y1) = ln(P_i^0) + ln(x_i*gamma_i),
%     gamma_i cannot be negative during optimizations
%     102 uses MultiStart() for mode 2
%   x0: initial guess for z_i, lnP_i
%   Q is a P x R matrix and contains loadings predicted by the IAST theory.
%
%   See Bai et al., Langmuir 28 (2012) 15566 for reference to equations
%   and symbols.
%
% See also IAST_func, lsqnonlin, optimset, interp1, ppval

function [Q_predicted, x, err, lnP0, psi, exitflags] = IAST_solve(M, S, varargin)
if nargin < 2 || rem(nargin,2) ~= 0
    error('IAST_solve:IncorrectNumberArguments', 'Number of arguments incorrect');
end

[ndata, N] = size(M);
if N < 2
    error('IAST_solve:TooFewComponents', 'System must have at least two components');
end

options_0 = optimset('FinDiffType', 'central', 'FunValCheck', 'on', 'MaxFunEvals', ndata*N*250, 'MaxIter', ndata*N*3, 'TolFun', 1e-6, 'TolX', 1e-6, 'Display', 'off');

names          = {'isotherm'; 'minlnP'; 'maxlnP'; 'EoS'; 'options'; 'mode'; 'x0'; 'x_lb'; 'x_ub'; 'tol'; 'EoS_deriv'; 'ads_pot'; 'inv_ads_pot'};
default_values = {        [];       [];       [];    []; options_0;      1;   [];     [];     [];  1e-5;          [];        []; []};
opt_args = process_variable_arguments(names, default_values, varargin);
isotherm = opt_args.('isotherm');
minlnP = opt_args.('minlnP');
maxlnP = opt_args.('maxlnP');
EoS = opt_args.('EoS');
options = opt_args.('options');
mode = opt_args.('mode');
x0 = opt_args.('x0');
x_lb = opt_args.('x_lb');
x_ub = opt_args.('x_ub');
tol = opt_args.('tol');
EoS_deriv =opt_args.('EoS_deriv');
ads_pot = opt_args.('ads_pot');
inv_ads_pot = opt_args.('inv_ads_pot');

if isempty(isotherm)
    [isotherm, minlnP] = fit_piecewise_polynomial(S);
end

if isempty(maxlnP)
    if isempty(S)
        error('IAST_solve:InsufficientParameters', 'Either S or maxlnP needs to be provided');
    end
    maxlnP = zeros(1, N);
    for i = 1:N
        maxlnP(i) = max(log(S{i}(:, 1)));
    end
end

if isempty(EoS)
    EoS = @(x)ones(1, N);
end

if isempty(x_lb)
    x_lb = zeros(1, N-1);
end

if isempty(x_ub)
    x_ub = ones(1, N-1);
end

options = optimset(options_0, options);

lnP_mixture = log(M);

if isempty(x0)
    z = ones(ndata, N-1) * 1/N;
    x0 = IAST_init(lnP_mixture, z, EoS([z,1e3*ones(ndata,1)]), mode);
end

if isempty(EoS_deriv)
    % central difference to calculate numerical derivative
    EoS_deriv = @(x)(sum((log(EoS([x(1:end-1), x(end)+tol]))-log(EoS([x(1:end-1), x(end)-tol])))/2/tol.*x(1:end-1)));
end

Q_tot = zeros(ndata, 1);
Z = zeros(ndata, N);
lnP0 = zeros(ndata, N);
psi = zeros(ndata, N);
exitflags = zeros(ndata, 1);
if mode == 1 || mode == 2 || mode == 102
    err = zeros(ndata, 2*N-1);
    x = zeros(ndata, 2*N-1);
elseif mode == -1 || mode == -2
    err = zeros(ndata, N);
    x = zeros(ndata, N);
else
    error('IAST_solve:UnknownMode','Mode parameter outside known choices.')
end

for i = 1 : ndata  % mixture partial pressures
    if ~strcmp(optimget(options, 'Display'), 'off')
        fprintf('\n======Data point %d ======\n', i);
    end
    func = @(x)IAST_func(x, lnP_mixture(i, :), 'isotherm', isotherm, 'minlnP', minlnP, 'ads_pot', ads_pot, 'inv_ads_pot', inv_ads_pot, 'EoS', EoS, 'mode', mode, 'tol', tol);
    if mode == -1 || mode == -2
        lb = [x_lb, 0];
        ub = [x_ub, Inf];
    elseif mode == 1 || mode == 2 || mode == 102
        lb = [x_lb, minlnP];
        ub = [x_ub, 2*maxlnP];
    end
    if mode == 102
        prob = createOptimProblem('lsqnonlin', 'objective', func, 'x0', x0(i, :), ...
            'lb', lb, 'ub', ub, 'options', options);
        ms = MultiStart('Display', 'iter', 'TolFun', 1e-3, 'TolX', 1e-4, 'UseParallel', 'always');
        [x1, fval, exitflags(i), output, allmins] = run(ms, prob, 10);
    else
        [x1, resnorm, residual, exitflags(i), output, lambda, jacobian] = lsqnonlin(func, x0(i, :), lb, ub, options);
    end

    x(i, :) = x1;
    [err(i, :), lnP0(i, :), psi(i, :)] = func(x1);
    if strcmp(optimget(options, 'Display'), 'final')
        disp('err(i, :): ')
        disp(err(i, :));
    end
    Z(i, 1:N-1) = x1(1:N-1);
    Z(i, N) = 1 - sum(Z(i, 1:N-1));
    for j = 1 : N
        Q_tot(i) = Q_tot(i) + Z(i, j) / isotherm{j}(lnP0(i, j));
    end
    Q_tot(i) = 1 / (Q_tot(i) + EoS_deriv([Z(i, 1:N-1), psi(i, N)]));
end

Q_predicted = Q_tot .* Z;

end