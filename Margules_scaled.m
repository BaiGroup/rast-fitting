function gamma = Margules_scaled(coeff, x)
% Margules acitivity model with scaled \ln \gamma_{\text{water}}
% coeff(1): A21, coeff(2): A12, coeff(3): scaling factor
% x(j, i), 1<=i<=N-1: mole fraction of component i for j-th data point

gamma = Margules(coeff(1:2), x);
gamma = [exp(coeff(3)*log(gamma(:,1))), gamma(:,2)];

end