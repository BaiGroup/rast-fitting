function [err, Q_predicted, lnP0, psi] = RAST_obj(x, isotherm, Q, EoS_deriv, varargin)
% Real Adsorbed Solution Theory objective function
% x((i-1)*N+j), 1<=j<=N: ln p_j^0 for i-th data point
% x(ndata*N+(i-1)*(N-1)+j), 1<=j<=N-1: z_j for component j in the absorbed phase
% x(ndata*(2N-1)+j): j-th coefficient for EoS
% isotherm(i), 1<=i<=N: function handle for isotherm i, Q(lnP)

if nargin < 2 || rem(nargin,2) ~= 0
    error('RAST_obj:Number of arguments incorrect');
end

names          = {'minlnP'; 'tol'; 'ads_pot'};
default_values = {      [];  1e-5; []};
opt_args = process_variable_arguments(names, default_values, varargin);
minlnP = opt_args.('minlnP');
tol = opt_args.('tol');
ads_pot = opt_args.('ads_pot');

if isempty(ads_pot)
    ads_pot = cell(1, N);
end

[ndata, N] = size(Q);

coeff = x(ndata*(2*N-1)+1:end);
lnP0 = reshape(x(1:ndata*N), N, ndata)';
psi = zeros(ndata, N);
z = zeros(ndata, N);
z(:, 1:N-1) = reshape(x(ndata*N+1:ndata*(2*N-1)), N-1, ndata)';
z(:, N) = 1 - sum(z(:, 1:N-1), 2);

Q_tot = zeros(ndata, 1);
for i = 1 : ndata
    for j = 1 : N
        if isempty(ads_pot{j})
            if isempty(isotherm{j}) || isempty(minlnP)
                error('RAST_obj:isotherm and minlnP must be provided without ads_pot');
            end
            ads_pot{j} = @(y)adsorption_potential(y, isotherm{j}, minlnP(j), tol);
        end
        psi(i, j) = ads_pot{j}(lnP0(i, j));
        Q_tot(i) = Q_tot(i) + z(i, j) / isotherm{j}(lnP0(i, j));
    end
    Q_tot(i) = 1 / (Q_tot(i) + EoS_deriv(coeff, [z(i, 1:N-1), psi(i, N)]));
end
Q_predicted = Q_tot .* z;
err = sum(((Q_predicted - Q)./Q).^2, 'all');

end