function err=IAST_func(x,isotherm,minlnP,lnP,gamma,tol,mode)
% AST equations for N-component systems
% x(i), 1<=i<=N: ln p_i^0
% x(N+i), 1<=i<=N-1: z_i for component i in the absorbed phase
% isotherm(i), 1<=i<=N: function handle for isotherm i, Q(lnP)
% minlnP(i), 1<=i<=N: lnP at which Q is zero
% lnP(i), 1<=i<=N: bulk phase pressure
% gamma(i): activity coefficient for component i; default is 1
% tol: precision
% mode: 1 uses P*y1/P_i^0 - x_i*gamma_i = 0
%       2 uses ln(P*y1) = ln(P_i^0) + ln(x_i*gamma_i), gamma_i cannot be
%       negative during optimizations

N=length(minlnP);
z=[x(N+1:2*N-1), 1-sum(x(N+1:2*N-1))];

if nargin<5 || isempty(gamma)
    gamma=ones(1,N);
end

if nargin < 6 || isempty(tol)
    tol = 1e-5;
end

if nargin < 7
    mode = 1;
end

err=zeros(1,2*N-1);
funcN=isotherm{N};
if x(N) < minlnP(N) + tol
    tN = 0;
else
    tN = trapz(minlnP(N):tol:x(N), funcN(minlnP(N):tol:x(N)));
end
for i=1:N-1
    funcI=isotherm{i};
    if x(i) < minlnP(i) + tol
        t = 0;
    else
        t = trapz(minlnP(i):tol:x(i), funcI(minlnP(i):tol:x(i)));
    end
    err(i) = tN - t;
end

if mode == 1
    err(N:2*N-1)=exp(lnP-x(1:N))-z.*gamma;
else
    err(N:2*N-1)=lnP-x(1:N)-log(z.*gamma);
end
end