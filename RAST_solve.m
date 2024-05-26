function [predicted, x, err, lnP0, psi, Q_IAST, x_IAST, err_IAST, lnP0_IAST, psi_IAST] = RAST_solve(M, S, varargin)
% IAST solver for N-component systems
% M{j}(i, 1:2): component j, i-th partial pressure & loading
% S{j}(i, 1:2): component j, i-th pressure & loading
% EoS: function handle that computes the activity coefficients
%      [\gamma_1, ..., \gamma_N] = EoS([z_1, z_2, ..., z_N])
% options: argument for fmincon()
% mode: 1 or 2 finds the best-fit activity coefficients for each data point
%       1 uses P*y1/P_i^0 - x_i*gamma_i = 0
%       2 uses ln(P*y1) = ln(P_i^0) + ln(x_i*gamma_i), gamma_i cannot be
%       negative during optimizations
%       >2 solves for the EoS parameters that best reproduces the entire
%       binary isotherm
%       3 Minimize the squared errors in loadings with the AST
%       equations as constraints
%       4 Minimize the squared errors in AST equations AND loadings
%       5 Minimize the squared errors in loadings with AST equations
%       solved by IAST_solve() at each iteration. Computationally extremely
%       expensive.
%       6 Minimize the squared errors in ln(P) and S_1j
% x0: initial guess for EoS parameters for mode > 1
% N_EoS_param: number of parameters in the activity model

if nargin < 2 || rem(nargin,2) ~= 0
    error('RAST_solve:Number of arguments incorrect');
end

N = length(M);
ndata = size(M{1}, 1);

options_0 = optimset('FinDiffType', 'central', 'FunValCheck', 'on', 'MaxFunEvals', 6000, 'MaxIter', 100, 'TolFun', 1e-6, 'TolX', 1e-6, 'Display', 'iter');
names          = {'isotherm'; 'minlnP'; 'maxlnP'; 'EoS'; 'options'; 'mode'; 'x0'; 'x_lb'; 'x_ub'; 'tol'; 'N_EoS_param'; 'EoS_deriv'; 'ads_pot'; 'inv_ads_pot'; 'C_ub'; 'g_lb'; 'g_ub'; 'batch'; 'skip'};
default_values = {        [];       [];       [];    []; options_0;      1;   [];      0;      1;  1e-5;            [];          [];        [];            [];    Inf;   -Inf;    Inf; 2*ndata; ndata};
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
N_EoS_param = opt_args.('N_EoS_param');
EoS_deriv = opt_args.('EoS_deriv');
ads_pot = opt_args.('ads_pot');
inv_ads_pot = opt_args.('inv_ads_pot');
C_ub = opt_args.('C_ub');
g_lb = opt_args.('g_lb');
g_ub = opt_args.('g_ub');
batch = opt_args.('batch');
skip = opt_args.('skip');

options = optimset(options_0, options);

if isempty(x0)
    noX0 = true;
    if mode > 2 && mode < 100 && isempty(N_EoS_param)
        error('RAST_solve:NotEnoughParameters', 'N_EoS_param must be provided when x0 is empty')
    end
else
    noX0 = false;
end

if isempty(x_lb)
    x_lb = zeros(1, N-1);
end

if isempty(x_ub)
    x_ub = ones(1, N-1);
end

lnP = zeros(ndata, N);
Q = zeros(ndata, N);
for i = 1 : N  % components
    lnP(:, i) = log(M{i}(:, 1));
    Q(:, i) = M{i}(:, 2);
end

if isempty(isotherm)
    [isotherm, minlnP, maxlnP, ads_pot, inv_ads_pot] = fit_piecewise_polynomial(S);
end

if isempty(maxlnP)
    if isempty(S)
        error('RAST_solve:InsufficientParameters', 'Either S or maxlnP needs to be provided');
    end
    maxlnP = zeros(1, N);
    for i = 1:N
        maxlnP(i) = max(log(S{i}(:, 1)));
    end
end

if isempty(EoS_deriv)
    % central difference to calculate numerical derivative
    EoS_deriv = @(y)(sum((log(EoS(coeff, [y(1:end-1), y(end)+tol]))-log(EoS(coeff, [y(1:end-1), y(end)-tol])))/2/tol.*y(1:end-1)));
end

if noX0 && mode ~= 5 && mode ~= 6
    % Initialize using IAST solutions and activity coefficients of one
    % (probably by setting EoS parameters to zero)
    if ~strcmp(optimget(options, 'Display'), 'off')
        disp('Solve IAST for initial guess')
        options_IAST = optimset('Display', 'iter');
        [Q_IAST, x_IAST, err_IAST, lnP0_IAST, psi_IAST] = IAST_solve(exp(lnP), [], 'isotherm', isotherm, 'minlnP', minlnP, 'maxlnP', maxlnP, 'ads_pot', ads_pot, 'inv_ads_pot', inv_ads_pot, 'options', options_IAST, 'mode', 1);
    end
else
    Q_IAST = [];
    x_IAST = [];
    err_IAST = [];
    lnP0_IAST = [];
    psi_IAST = [];
end

lnP0 = zeros(ndata, N);
psi = zeros(ndata, N);

if mode == 3
    if noX0
        x0 = [reshape(x_IAST(:, N:end)', 1, []), reshape(x_IAST(:, 1:N-1)', 1, []), zeros(1, N_EoS_param)];
    else
        N_EoS_param = length(x0) - (2*N-1)*ndata;
    end
    objFunc = @(x)RAST_obj(x, isotherm, Q, EoS_deriv, 'minlnP', minlnP, 'ads_pot', ads_pot, 'tol', tol);
    conFunc = @(x)RAST_const(x, isotherm, lnP, EoS, 'minlnP', minlnP, 'ads_pot', ads_pot, 'inv_ads_pot', inv_ads_pot, 'tol', tol);
    lb = [repmat(minlnP, 1, ndata), repmat(x_lb, 1, ndata), g_lb*ones(1, N_EoS_param-1), 0];
    ub = [repmat(maxlnP, 1, ndata), repmat(x_ub, 1, ndata), g_ub*ones(1, N_EoS_param-1), C_ub];
    [x, fval, exitflag, output, lambda, grad, hessian] = fmincon(objFunc, x0, [], [], [], [], lb, ub, conFunc, options);
    [err, predicted, lnP0, psi] = objFunc(x);
elseif mode == 4
    if noX0
        x0 = [reshape(x_IAST(:, N:2*N-1)', 1, []), zeros(1, N_EoS_param)];
    else
        N_EoS_param = length(x0)-N*ndata;
    end
    x = x0;
    options_b = optimset(options, 'TolFun', 1e-4, 'TolX', 1e-3, 'MaxIter', 50);  %'OutputFcn', @outfun);
    for i = 1:skip:ndata-batch+1
        for j = 1:N
            Mb{j}=M{j}(i:i+batch-1, :);
        end
        x0b = [x(N*(i-1)+1:N*(i-1+batch)), x(end-N_EoS_param+1:end)]
        objFunc = @(y)RAST_func(y, isotherm, Mb, EoS, 'minlnP', minlnP, 'ads_pot', ads_pot, 'inv_ads_pot', inv_ads_pot, 'EoS_deriv', EoS_deriv, 'tol', tol);
        lb = [repmat(minlnP,1,batch), -Inf*ones(1,N_EoS_param-1), 0];
        [xb, fval, exitflag, output, lambda, grad, hessian] = fmincon(objFunc, x0b, [], [], [], [], lb, [], [], options_b);
        x(N*(i-1)+1:N*(i-1+batch)) = xb(1:end-N_EoS_param);
        x(end-N_EoS_param+1:end) = (x(end-N_EoS_param+1:end)*(i-1) + xb(end-N_EoS_param+1:end))/i;
    end
    objFunc = @(y)RAST_func(y, isotherm, M, EoS, 'minlnP', minlnP, 'ads_pot', ads_pot, 'inv_ads_pot', inv_ads_pot, 'EoS_deriv', EoS_deriv, 'tol', tol);
    lb = [repmat(minlnP, 1, ndata), g_lb*ones(1, N_EoS_param-1), 0];
    ub = [repmat(2*maxlnP, 1, ndata), g_ub*ones(1, N_EoS_param-1), C_ub];
    [x, fval, exitflag, output, lambda, grad, hessian] = fmincon(objFunc, x, [], [], [], [], lb, ub, [], options);
    [err, predicted, lnP0, psi] = objFunc(x);
elseif mode == 5
    if noX0
        x0 = zeros(1, N_EoS_param);  % many EoS return 1 as activity coefficients when all parameters are zero
    else
        N_EoS_param = length(x0);
    end
    objFunc = @(x)RAST_func_IAST_solve(x, isotherm, M, EoS, 'minlnP', minlnP, 'maxlnP', maxlnP, 'ads_pot', ads_pot, 'inv_ads_pot', inv_ads_pot, 'EoS_deriv', EoS_deriv, 'tol', tol);
    lb = [g_lb*ones(1, N_EoS_param-1), 0];
    ub = [g_ub*ones(1, N_EoS_param-1), C_ub];
    [x, fval, exitflag, output, lambda, grad, hessian] = fmincon(objFunc, x0, [], [], [], [], lb, ub, [], options);
    [err, Q_IAST, x_IAST, err_IAST, lnP0_IAST, psi_IAST] = objFunc(x);
    predicted = Q_IAST;
    lnP0 = lnP0_IAST;
    psi = psi_IAST;
elseif mode == 6
    if noX0
        x0 = zeros(1, N_EoS_param);  % many EoS return 1 as activity coefficients when all parameters are zero
    else
        N_EoS_param = length(x0);
    end
    objFunc = @(x)RAST_solve_for_fugacity(x, M, EoS, 'isotherm', isotherm, 'minlnP', minlnP, 'ads_pot', ads_pot, 'inv_ads_pot', inv_ads_pot, 'EoS_deriv', EoS_deriv, 'tol', tol);
    lb = [g_lb*ones(1, N_EoS_param-1), 0];
    ub = [g_ub*ones(1, N_EoS_param-1), C_ub];
    [x, fval, exitflag, output, lambda, grad, hessian] = fmincon(objFunc, x0, [], [], [], [], lb, ub, [], options);
    [err, predicted, lnP0, psi] = objFunc(x);
elseif mode == 1 || mode == 2 || mode == 102
    x = zeros(ndata, 2*N);
    predicted = zeros(ndata, N);
    psi = zeros(ndata, N);
    err = zeros(ndata, 3*N-1);
    if noX0
        x0 = [x_IAST(:, N:2*N-1), ones(ndata, N)];
    end
    for i = 1 : ndata  % mixture partial pressures
        fprintf('\n======Data point %d ======\n', i);
        func = @(x)RAST_func_per_point(x, isotherm, lnP(i, :), Q(i, :), 'minlnP', minlnP, 'ads_pot', ads_pot, 'inv_ads_pot', inv_ads_pot, 'mode', mode, 'tol', tol);
        lb = [minlnP, g_lb*ones(1, N)];
        ub = [2*maxlnP, g_ub*ones(1, N)];
        [x(i, :), resnorm, residual, exitflag, output, lambda, jacobian] = lsqnonlin(func, x0(i, :), lb, ub, options);
        [err(i, :), predicted(i, :), lnP0(i, :), psi(i, :)] = func(x(i,:));
    end
else
    error('RAST_solve:UnknownMode','Mode parameter outside known choices.')
end

    function stop = outfun(x, optimValues, state)
        stop = false;
        disp(state);
        disp(x);
    end

end