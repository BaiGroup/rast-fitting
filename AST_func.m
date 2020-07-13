function y=AST_func(x,isotherm,minlnP,lnP,tol,EOS)
% AST equations for N-component systems
% x(i), 1<=i<=N: ln p_i^0
% x(i), N+1<=i<=2N-1: z_{i-N} for component i-N the absorbed phase
% isotherm(i), 1<=i<=N: function handle for isotherm i
% minlnP(i), 1<=i<=N: lower limit of integration,
%       i.e., pressure at which the amount of adsorption is zero
% lnP(i), 1<=i<=N: bulk phase pressure
% tol: precision

N=length(minlnP);
xN=1-sum(x(N+1:end));
if nargin<6
    gamma=ones(N,1);
else
    gamma=EOS([x(N+1:end);xN]');
end
y=zeros(2*N-1,1);
funcN=isotherm{N};
for i=1:N-1
    funcI=isotherm{i};
    y(i)=trapz(minlnP(N):tol:x(N),funcN(minlnP(N):tol:x(N)))-trapz(minlnP(i):tol:x(i),funcI(minlnP(i):tol:x(i)));
    y(i+N-1)=exp(lnP(i)-x(i))-x(i+N)*gamma(i);  %% P*y1/P_i^0 - x_i*gamma_i normalizes residual
end
y(2*N-1)=exp(lnP(N)-x(N))-xN*gamma(N);

end