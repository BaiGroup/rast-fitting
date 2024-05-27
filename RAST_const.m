function [c, ceq] = RAST_const(x, isotherm, lnP, EoS, varargin)
% Real Adsorbed Solution Theory constraint function
% x((i-1)*N+j), 1<=j<=N: ln p_j^0 for i-th data point
% x(ndata*N+(i-1)*(N-1)+j), 1<=j<=N-1: z_j for component j in the absorbed phase
% x(ndata*(2N-1)+j): j-th coefficient for EoS
% isotherm(i), 1<=i<=N: function handle for isotherm i, Q(lnP)
% minlnP(i), 1<=i<=N: lnP at which Q is zero
% EoS: function handle that computes the activity coefficients
%      [\gamma_1, ..., \gamma_N] = EoS([coeff_1, ..., coeff_M], [z_1, ..., z_N-1, \Psi])

if nargin < 4 || rem(nargin,2) ~= 0
    error('RAST_const:Number of arguments incorrect');
end

names          = {'minlnP'; 'ads_pot'; 'inv_ads_pot'; 'tol'};
default_values = {      [];        [];            []; 1e-5};
opt_args = process_variable_arguments(names, default_values, varargin);
minlnP = opt_args.('minlnP');
ads_pot = opt_args.('ads_pot');
inv_ads_pot = opt_args.('inv_ads_pot');
tol = opt_args.('tol');

[ndata, N] = size(lnP);

coeff = x(ndata*(2*N-1)+1:end);
lnP0 = reshape(x(1:ndata*N), N, ndata)';
z = reshape(x(ndata*N+1:ndata*(2*N-1)), N-1, ndata)';

EoS_with_coeff = @(y)EoS(coeff, y);

err = zeros(ndata, 2*N-1);
for i = 1 : ndata
    err(i, :) = IAST_func([z(i, :), lnP0(i, :)], lnP(i, :), 'mode', 1, 'isotherm', isotherm, 'minlnP', minlnP, 'tol', tol, 'EoS', EoS_with_coeff, 'ads_pot', ads_pot, 'inv_ads_pot', inv_ads_pot);
end

c = [];
ceq = reshape(err, [], 1);

end