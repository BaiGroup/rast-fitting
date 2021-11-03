function gamma=Margules_ads(coeff,x)
% Margules acitivity model for adsorption
% coeff(1): A21, coeff(2): A12, coeff(3): C
% x(j, i), 1<=i<=N-1: mole fraction of component i for j-th data point
% x(j, N): adsorption potential

[ndata,N]=size(x);
if N ~= 2
    error('Margules_ads:DimensionError','Only works for binary mixture');
end
z=[x(:, 1:end-1),1-sum(x(:, 1:end-1),2)];

fac = 1 - exp(-coeff(3)*x(:, N));

gamma=[exp((z(:, 2).^2.*(coeff(2)+2*(coeff(1)-coeff(2))*z(:, 1))).*fac),exp((z(:, 1).^2.*(coeff(1)+2*(coeff(2)-coeff(1))*z(:, 2))).*fac)];
% gamma(isinf(gamma)) = 1e3;

end