%% MFI-1/298K data
load('MFI-1-298K.mat');

% Calculate vapor pressure for high pressure liquid water
p_sat_H2O=exp(1.50116)*1e3;  % T=298K
V_H2O=18.015e-6/0.991665;  % g/mL -> m^3/mol
water(:,1)=exp((water(:,1)-p_sat_H2O)*V_H2O/8.314/298)*p_sat_H2O;

%% MFI-1/323K data
load('MFI-1-323K.mat');

% Calculate vapor pressure for high pressure liquid water
p_sat_H2O=18.165e3;  % T=323K
V_H2O=18.46e-6;  % m^3/mol
% water(:,1)=exp((water(:,1)-p_sat_H2O)*V_H2O/8.314/323)*p_sat_H2O;

%% water/methanol
options=optimset('Display','iter');
S={[water(:,1),water(:,2)],[methanol(:,1),methanol(:,2)]};
M=[methanol_water(:,1),methanol_water(:,2)];
[Q, x, err, lnP0, psi] = IAST_solve(M, S, 'mode', 1, 'options', options);

semilogx(M(:,2),methanol_water(:,4),'rs',M(:,2),methanol_water(:,3),'bo',M(:,2),Q(:,2),'md',M(:,2),Q(:,1),'c^');
xlabel('$p$ [Pa]'); ylabel('$N$ [mol/kg]');
legend('sim\_MeOH','sim\_H2O','IAST\_MeOH','IAST\_H2O','Location','NorthEastOutside');

%% water/ethanol
options=optimset('Display','iter');
S={[water(:,1),water(:,2)],[ethanol(:,1),ethanol(:,2)]};
M=[ethanol_water(:,1),ethanol_water(:,2)];
[Q, x, err, lnP0, psi] = IAST_solve(M, S, 'mode', 1, 'options', options);
% [isotherm, minlnP, maxlnP, ads_pot, inv_ads_pot] = fit_piecewise_polynomial(S(1));
% [isotherm(2), minlnP(2), maxlnP(2), ads_pot(2), inv_ads_pot(2)] = fit_piecewise_polynomial(S(2));
% [Q, x, err, lnP0, psi] = IAST_solve(M, [], 'mode', 1, 'tol', 1e-6, 'isotherm', isotherm, 'minlnP', minlnP, 'maxlnP', maxlnP, 'ads_pot', ads_pot, 'inv_ads_pot', inv_ads_pot, 'EoS_deriv', @(x,y)0, 'options', options);

%% Plotting
semilogx(M(:,2),ethanol_water(:,4),'rs',M(:,2),ethanol_water(:,3),'bo',M(:,2),Q(:,2),'md',M(:,2),Q(:,1),'c^');
xlabel('$p$ [Pa]'); ylabel('$N$ [mol/kg]');
legend('sim\_EtOH','sim\_H2O','IAST\_EtOH','IAST\_H2O','Location','NorthEastOutside');

%% water/methanol/ethanol
options=optimset('Display','iter');
S={[water(:,1),water(:,2)],[methanol(:,1),methanol(:,2)],[ethanol(:,1),ethanol(:,2)]};
M=[ternary(:,1),ternary(:,2),ternary(:,3)];
[Q, x, err, lnP0, psi] = IAST_solve(M, S, 'mode', 1, 'options', options);

plot([0,1],[0,1],'k');hold on
plot(ternary(:,5)/(ternary(:,4)+ternary(:,5)+ternary(:,6)),Q(:,2)/(Q(:,1)+Q(:,2)+Q(:,3)),'mv');hold on
plot(ternary(:,6)/(ternary(:,4)+ternary(:,5)+ternary(:,6)),Q(:,3)/(Q(:,1)+Q(:,2)+Q(:,3)),'ch');
xlabel('$x_{\textrm{sim}}$'); ylabel('$x_{\textrm{IAST}}$');