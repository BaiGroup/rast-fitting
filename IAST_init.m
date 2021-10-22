function x0 = IAST_init(lnP, z, gamma)
[ndata, N] = size(lnP);

if nargin<3 || isempty(gamma)
    gamma=ones(ndata,N);
end

x0 = zeros(ndata, 2*N);

if nargin < 2
    x0(1:ndata, N+1:2*N) = 1/N;  % equimolar mixtures
else
    x0(1:ndata, N+1:2*N) = [z, 1-sum(z,2)];
end

x0(1:ndata, 1:N) = lnP - log(x0(1:ndata, N+1:2*N)) - log(gamma);  % lnP_i^0 = lnP_i - lnz_i -ln(\gamma_i)

x0 = x0(1:ndata, 1:end-1);
end