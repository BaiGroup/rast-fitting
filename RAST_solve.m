function [Q_predicted, x, Q_IAST, x_IAST, F_IAST] = RAST_solve(M, S, EoS, options, mode, x0, N_EoS_param)
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

options_0 = optimset('FinDiffType', 'central', 'FunValCheck', 'on', 'MaxFunEvals', 6000, 'MaxIter', 100, 'TolFun', 1e-6, 'TolX', 1e-6, 'Display', 'iter');
if nargin < 4 || isempty(options)
    options = options_0;
else
    options = optimset(options_0, options);
end

if nargin < 5 || isempty(mode)
    mode = 1;
end

if nargin < 6 || isempty(x0)
    noX0 = true;
    if mode > 2 && (nargin < 7 || isempty(N_EoS_param))
        error('RAST_solve:NotEnoughParameters', 'N_EoS_param must be provided when x0 is empty')
    end
else
    noX0 = false;
end

N = length(M);
ndata = size(M{1}, 1);
lnP = zeros(ndata, N);
Q = zeros(ndata, N);
for i = 1 : N  % components
    lnP(:, i) = log(M{i}(:, 1));
    Q(:, i) = M{i}(:, 2);
end

[isotherm, minlnP] = fit_isotherm(S);

if noX0 && mode ~= 5
    % Initialize using IAST solutions and activity coefficients of one
    % (probably by setting EoS parameters to zero)
    disp('Solve IAST for initial guess')
    [Q_IAST, x_IAST, F_IAST] = IAST_solve(exp(lnP), [], isotherm, minlnP, [], options, 1);
else
    Q_IAST = [];
    x_IAST = [];
    F_IAST = [];
end

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
        x0 = [reshape(x_IAST(:, 1:N)', 1, []), zeros(1, N_EoS_param)];
    end
    objFunc = @(x)RAST_func(x, isotherm, minlnP, M, EoS);
    lb = [repmat(minlnP,1,ndata), -Inf*ones(1,length(x0)-N*ndata)];
    [x,fval,exitflag,output,lambda,grad,hessian]=fmincon(objFunc,x0,[],[],[],[],lb,[],[],options);
    [err, Q_predicted] = objFunc(x);
elseif mode == 5
    if noX0
        x0 = zeros(1, N_EoS_param);  % many EoS return 1 as activity coefficients when all parameters are zero
    end
    objFunc=@(x)RAST_func_IAST_solve(x, isotherm, minlnP, M, EoS);
    [x,fval,exitflag,output]=fminsearch(objFunc,x0,options);
    [err, Q_predicted] = objFunc(x);
else
    x = zeros(ndata, 3*N-1);
    Q_predicted = zeros(ndata, N);
    if noX0
        x0 = [x_IAST, ones(ndata, N)];
    end
    for i = 1 : ndata  % mixture partial pressures
        i
        func = @(x)RAST_func_per_point(x, isotherm, minlnP, lnP(i, :), Q(i, :), mode);
        if mode == 1
            [x(i,:),fval,exitflag,output,jacobian] = fsolve(func, x0(i, :), options);
        else
            lb = [minlnP, zeros(1, N-1), -Inf*zeros(1, N)];
            ub = [Inf*ones(1, N), ones(1, N-1), Inf*zeros(1, N)];
            [x(i,:),resnorm,residual,exitflag,output,lambda,jacobian] = lsqnonlin(func, x0(i, :), lb, ub, options);
        end
        [err, Q_predicted(i, :)] = func(x(i,:));
    end
end

end