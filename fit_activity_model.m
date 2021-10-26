function coeff = fit_activity_model(EoS, N_EoS_param, z, gamma, options)
% Fit numerical interpolants for N single-component isotherms
% EoS: function handle that computes the activity coefficients
%      [\gamma_1, ..., \gamma_N] = EoS([coeff_1, ..., coeff_M], [z_1, ..., z_N])
% N_EoS_param: number of parameters in the activity model
% z(i), 1<=i<=N-1: z_i for component i in the absorbed phase
% gamma(i): activity coefficient for component i; default is 1
% options: argument for fminsearch()

options_0 = optimset('FinDiffType', 'central', 'FunValCheck', 'on', 'MaxFunEvals', 600, 'MaxIter', 100, 'Display', 'iter');
if nargin < 5 || isempty(options)
    options = options_0;
else
    options = optimset(options_0, options);
end

x0 = ones(1, N_EoS_param);
% func = @(x)(EoS(x, z)-gamma);
func = @(x)((EoS(x, z)-gamma)./gamma);
% func = @(x)(log(EoS(x, z)./gamma));
[coeff,fval,exitflag,output,jacobian] = fsolve(func, x0, options);

end