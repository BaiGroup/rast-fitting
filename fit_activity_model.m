function coeff = fit_activity_model(EoS, N_EoS_param, z, gamma, varargin)
% Fit numerical interpolants for N single-component isotherms
% EoS: function handle that computes the activity coefficients
%      [\gamma_1, ..., \gamma_N] = EoS([coeff_1, ..., coeff_M], [z_1, ..., z_N-1, \Psi])
% N_EoS_param: number of parameters in the activity model
% z(i), 1<=i<=N-1: z_i for component i in the absorbed phase
% gamma(i): activity coefficient for component i; default is 1
% options: argument for lsqnonlin()

if nargin < 4 || rem(nargin,2) ~= 0
    error('fit_activity_model:Number of arguments incorrect');
end

options_0 = optimset('FinDiffType', 'central', 'FunValCheck', 'on', 'MaxFunEvals', N_EoS_param*200, 'MaxIter', N_EoS_param*35, 'Display', 'iter');
names          = {'options'; 'C_ub'; 'g_lb'; 'g_ub'};
default_values = {options_0;    Inf;   -Inf;    Inf};
opt_args = process_variable_arguments(names, default_values, varargin);
options = opt_args.('options');
C_ub = opt_args.('C_ub');
g_lb = opt_args.('g_lb');
g_ub = opt_args.('g_ub');

options = optimset(options_0, options);
lb = [g_lb*ones(1, N_EoS_param-1), 0];
ub = [g_ub*ones(1, N_EoS_param-1), C_ub];

x0 = ones(1, N_EoS_param);
% func = @(x)(EoS(x, z)-gamma);
func = @(x)((EoS(x, z)-gamma)./gamma);
% func = @(x)(log(EoS(x, z)./gamma));
[coeff, resnorm, residual, exitflags, output, lambda, jacobian] = lsqnonlin(func, x0, lb, ub, options)

end