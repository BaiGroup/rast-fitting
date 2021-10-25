function [Q_predicted, x, err, lnP0, psi, Q_IAST, x_IAST, err_IAST, lnP0_IAST, psi_IAST] = RAST_solve(M, S, varargin)
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
%       3 Minimize the least-square error in loadings with the AST
%       equations as constraints
%       4 Minimize the least-square error in AST equations AND loadings
%       5 Minimize the least-square error in loadings with AST equations
%       solved by IAST_solve() at each iteration. Computationally extremely
%       expensive.
% x0: initial guess for EoS parameters for mode > 1
% N_EoS_param: number of parameters in the activity model

if nargin < 2 || rem(nargin,2) ~= 0
    error('RAST_solve:Number of arguments incorrect');
end

N = length(M);
ndata = size(M{1}, 1);

options_0 = optimset('FinDiffType', 'central', 'FunValCheck', 'on', 'MaxFunEvals', 6000, 'MaxIter', 100, 'TolFun', 1e-4, 'TolX', 1e-4, 'Display', 'iter');
names          = {'isotherm'; 'minlnP'; 'EoS'; 'options'; 'mode'; 'x0'; 'N_EoS_param'; 'EoS_deriv'; 'ads_pot'; 'inv_ads_pot'; 'batch'; 'skip'};
default_values = {        [];       [];    []; options_0;      1;   [];            [];          [];        [];            []; 2*ndata; ndata};
opt_args = process_variable_arguments(names, default_values, varargin);
isotherm = opt_args.('isotherm');
minlnP = opt_args.('minlnP');
EoS = opt_args.('EoS');
options = opt_args.('options');
mode = opt_args.('mode');
x0 = opt_args.('x0');
N_EoS_param = opt_args.('N_EoS_param');
EoS_deriv = opt_args.('EoS_deriv');
ads_pot = opt_args.('ads_pot');
inv_ads_pot = opt_args.('inv_ads_pot');
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

lnP = zeros(ndata, N);
Q = zeros(ndata, N);
for i = 1 : N  % components
    lnP(:, i) = log(M{i}(:, 1));
    Q(:, i) = M{i}(:, 2);
end

if isempty(isotherm) || isempty(minlnP)
    [isotherm, minlnP, maxlnP] = fit_isotherm(S);
end

if noX0 && mode ~= 5
    % Initialize using IAST solutions and activity coefficients of one
    % (probably by setting EoS parameters to zero)
    if ~strcmp(optimget(options, 'Display'), 'off')
        disp('Solve IAST for initial guess')
        options_IAST = optimset('Display', 'iter');
        [Q_IAST, x_IAST, err_IAST, lnP0_IAST, psi_IAST] = IAST_solve(exp(lnP), [], 'isotherm', isotherm, 'minlnP', minlnP, 'ads_pot', ads_pot, 'inv_ads_pot', inv_ads_pot, 'options', options_IAST, 'mode', 1);
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
        x0 = [reshape(x_IAST(:, 1:N)', 1, []), reshape(x_IAST(:, N+1:end)', 1, []), zeros(1, N_EoS_param)];
    end
    objFunc = @(x)RAST_obj(x, isotherm, M);
    conFunc = @(x)RAST_const(x, isotherm, minlnP, M, EoS);
    lb = [repmat(minlnP,1,ndata), zeros(1, (N-1)*ndata), -Inf*ones(1,length(x0)-(2*N-1)*ndata)];
    ub = [Inf*ones(1,N*ndata), ones(1, (N-1)*ndata), Inf*ones(1,length(x0)-(2*N-1)*ndata)];
    [x,fval,exitflag,output,lambda,grad,hessian]=fmincon(objFunc,x0,[],[],[],[],lb,ub,conFunc,options);
    [err, Q_predicted] = objFunc(x);
elseif mode == 4
    if noX0
        x0 = [reshape(x_IAST(:, N:2*N-1)', 1, []), ones(1, N_EoS_param)];
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
        objFunc = @(y)RAST_func(y, isotherm, Mb, EoS, 'minlnP', minlnP, 'ads_pot', ads_pot, 'inv_ads_pot', inv_ads_pot, 'EoS_deriv', EoS_deriv);
        lb = [repmat(minlnP,1,batch), -Inf*ones(1,N_EoS_param-1), 0];
        [xb,fval,exitflag,output,lambda,grad,hessian]=fmincon(objFunc,x0b,[],[],[],[],lb,[],[],options_b);
        x(N*(i-1)+1:N*(i-1+batch)) = xb(1:end-N_EoS_param);
        x(end-N_EoS_param+1:end) = (x(end-N_EoS_param+1:end)*(i-1) + xb(end-N_EoS_param+1:end))/i
    end
    objFunc = @(y)RAST_func(y, isotherm, M, EoS, 'minlnP', minlnP, 'ads_pot', ads_pot, 'inv_ads_pot', inv_ads_pot, 'EoS_deriv', EoS_deriv);
    lb = [repmat(minlnP,1,ndata), -Inf*ones(1,N_EoS_param-1), 0];
    ub = [Inf*ones(1, ndata*N), Inf*ones(1,N_EoS_param-1), 1e-2];
    [x,fval,exitflag,output,lambda,grad,hessian]=fmincon(objFunc,x,[],[],[],[],lb,ub,[],options);
    [err, Q_predicted, lnP0, psi] = objFunc(x);
elseif mode == 5
    if noX0
        x0 = ones(1, N_EoS_param);  % many EoS return 1 as activity coefficients when all parameters are zero
    end
    objFunc = @(x)RAST_func_IAST_solve(x, isotherm, M, EoS, 'minlnP', minlnP, 'ads_pot', ads_pot, 'inv_ads_pot', inv_ads_pot, 'EoS_deriv', EoS_deriv);
%     [x, fval, exitflag, output] = fminsearch(objFunc, x0, options);
    lb = [-Inf*ones(1,N_EoS_param-1), 0];
    ub = [Inf*ones(1,N_EoS_param-1), 1e-2];
    [x,fval,exitflag,output,lambda,grad,hessian]=fmincon(objFunc,x0,[],[],[],[],lb,ub,[],options);
    [err, Q_predicted, x_IAST, err_IAST, lnP0_IAST, psi_IAST] = objFunc(x);
    lnP0 = lnP0_IAST;
    psi = psi_IAST;
elseif mode == 1 || mode == 2 || mode == 102
    x = zeros(ndata, 2*N);
    Q_predicted = zeros(ndata, N);
    psi = zeros(ndata, N);
    err = zeros(ndata, 3*N-1);
    if noX0
        x0 = [x_IAST(:, N:2*N-1), ones(ndata, N)];
    end
    for i = 1 : ndata  % mixture partial pressures
        fprintf('\n======Data point %d ======\n', i);
        func = @(x)RAST_func_per_point(x, isotherm, lnP(i, :), Q(i, :), 'minlnP', minlnP, 'ads_pot', ads_pot, 'inv_ads_pot', inv_ads_pot, 'mode', mode);
        if mode == 1
            [x(i,:),fval,exitflag,output,jacobian] = fsolve(func, x0(i, :), options);
        elseif mode == 2 || mode == 102
            lb = [minlnP, -Inf*zeros(1, N)];
            ub = [Inf*ones(1, N), Inf*zeros(1, N)];
            [x(i,:),resnorm,residual,exitflag,output,lambda,jacobian] = lsqnonlin(func, x0(i, :), lb, ub, options);
        end
        [err(i, :), Q_predicted(i, :), lnP0(i, :), psi(i, :)] = func(x(i,:));
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