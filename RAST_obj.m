function [err, Q_predicted] = RAST_obj(x, isotherm, M)
% Real Adsorbed Solution Theory objective function
% x((i-1)*N+j), 1<=j<=N: ln p_j^0 for i-th data point
% x(ndata*N+(i-1)*(N-1)+j), 1<=j<=N-1: z_j for component j in the absorbed phase
% x(ndata*(2N-1)+j): j-th coefficient for EoS
% isotherm(i), 1<=i<=N: function handle for isotherm i, Q(lnP)
% M{j}(i, 1:2): component j, i-th partial pressure & loading

N = length(M);
ndata = size(M{1}, 1);
lnP = zeros(ndata, N);
Q = zeros(ndata, N);
for i = 1 : N  % components
    lnP(:, i) = log(M{i}(:, 1));
    Q(:, i) = M{i}(:, 2);
end

lnP0 = reshape(x(1:ndata*N), N, ndata)';
z = zeros(ndata, N);
z(:, 1:N-1) = reshape(x(ndata*N+1:ndata*(2*N-1)), N-1, ndata)';
z(:, N) = 1 - sum(z(:, 1:N-1), 2);
Q_tot = zeros(ndata, 1);
for i = 1 : ndata
    for j=1:N
        Q_tot(i) = Q_tot(i) + z(i, j) / isotherm{j}(lnP0(i, j));
    end
    Q_tot(i) = 1 / Q_tot(i);
end
Q_predicted = Q_tot .* z;
err = sum(((Q_predicted - Q)./Q).^2, 'all');

end