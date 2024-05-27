function gamma=Margules(coeff,x)
% Margules acitivity model
% coeff(1): A21, coeff(2): A12
% x(j, i), 1<=i<=N-1: mole fraction of component i for j-th data point

[ndata,N]=size(x);
if N ~= 2
    error('Margules_ads:DimensionError','Only works for binary mixture');
end
z=[x(:, 1:end-1),1-sum(x(:, 1:end-1),2)];

% gamma=[exp(z(:, 2).^2.*(coeff(2)+2*(coeff(1)-coeff(2))*z(:, 1))),ones(ndata,1)];
gamma=[exp(z(:, 2).^2.*(coeff(2)+2*(coeff(1)-coeff(2))*z(:, 1))),exp(z(:, 1).^2.*(coeff(1)+2*(coeff(2)-coeff(1))*z(:, 2)))];

end