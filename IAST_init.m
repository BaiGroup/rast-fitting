function x0 = IAST_init(lnP, z, gamma, mode)
[ndata, N] = size(lnP);

if nargin<3 || isempty(gamma)
    gamma=ones(ndata,N);
end

if nargin < 4 || isempty(mode)
    mode = 1;
end

if nargin < 2 || isempty(z)
    x = 1/N*ones(ndata, N);  % equimolar mixtures
else
    x = [z, 1-sum(z,2)];
end

if mode == -1 || mode == -2
    x0 = [x(1:ndata,1:end-1), 50*ones(ndata,1)];  % Psi = 50
elseif mode == 1 || mode == 2 || mode == 102
    x0 = [x(1:ndata,1:end-1), lnP - log(x) - log(gamma)];  % lnP_i^0 = lnP_i - lnz_i -ln(\gamma_i)
else
    error('IAST_init:UnknownMode','Mode parameter outside known choices.')
end
end