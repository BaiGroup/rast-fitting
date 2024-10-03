function gamma=Margules_3C(coeff,x)
% Margules acitivity model for adsorption
% coeff(1): A21, coeff(2): A12, coeff(3): A31, coeff(4): A13,
% coeff(5): B21, coeff(6): B12, coeff(7): B32, coeff(8): B23,
% coeff(9): C31, coeff(10): C13, coeff(11): C32, coeff(12): C23,
% coeff(13): C
% x(j, i), 1<=i<=N-1: mole fraction of component i for j-th data point
% x(j, N): adsorption potential

[ndata,N]=size(x);
if N ~= 3
    error('Margules_ads:DimensionError','Only works for ternary mixture');
end
z=[x(:, 1:end-1),1-sum(x(:, 1:end-1),2)];

gamma=exp([z(:, 2).^2.*(coeff(2)+coeff(1)*z(:, 1)) + z(:, 3).^2.*(coeff(4)+coeff(3)*z(:, 1)), ...
    z(:, 1).^2.*(coeff(6)+coeff(5)*z(:, 2)) + z(:, 3).^2.*(coeff(8)+coeff(7)*z(:, 2)), ...
    z(:, 1).^2.*(coeff(10)+coeff(9)*z(:, 3)) + z(:, 2).^2.*(coeff(12)+coeff(11)*z(:, 3))]);

end