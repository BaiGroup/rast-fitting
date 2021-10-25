%% MFI-1/298K data
load('MFI-1-298K.mat');
set(0,'DefaultTextInterpreter','latex');
options=optimset('Display','iter');

% Calculate vapor pressure for high pressure liquid water
p_sat_H2O=exp(1.50116)*1e3;  % T=298K
V_H2O=18.015e-6/0.991665;  % g/mL -> m^3/mol
water(:,1)=exp((water(:,1)-p_sat_H2O)*V_H2O/8.314/298)*p_sat_H2O;

%% LTA-0/323K data
load('LTA-0-323K.mat');
set(0,'DefaultTextInterpreter','latex');
options=optimset('Display','iter');

% Calculate vapor pressure for high pressure liquid water
p_sat_H2O=18.165e3;  % T=323K
V_H2O=18.46e-6;  % m^3/mol
% water(:,1)=exp((water(:,1)-p_sat_H2O)*V_H2O/8.314/323)*p_sat_H2O;

%% FAU-2/323K data
load('FAU-2-323K.mat');
set(0,'DefaultTextInterpreter','latex');
options=optimset('Display','iter');

% Calculate vapor pressure for high pressure liquid water
p_sat_H2O=18.165e3;  % T=323K
V_H2O=18.46e-6;  % m^3/mol
% water(:,1)=exp((water(:,1)-p_sat_H2O)*V_H2O/8.314/323)*p_sat_H2O;

%% water/methanol
S={[water(:,1),water(:,2)],[methanol(:,1),methanol(:,2)]};
M=[methanol_water(:,6)*1000,methanol_water(:,3)*1000];
[Q, x, err, lnP0, psi] = IAST_solve(M, S, 'mode', 1, 'options', options);

semilogx(M(:,2),methanol_water(:,2),'rs',M(:,2),methanol_water(:,5),'bo',M(:,2),Q(:,2),'md',M(:,2),Q(:,1),'c^');
xlabel('$p$ [Pa]'); ylabel('$N$ [molec/uc]');
legend('sim\_MeOH','sim\_H2O','IAST\_MeOH','IAST\_H2O','Location','NorthEastOutside');

%% water/ethanol
S={[water(:,1),water(:,2)],[ethanol(:,1),ethanol(:,2)]};
M=[ethanol_water(:,6)*1000,ethanol_water(:,3)*1000];
[isotherm, ads_pot, inv_ads_pot] = fit_Langmuir_Sips(S(1));  % fit water to Langmuir-Sips isotherms
[isotherm2, minlnP2] = fit_isotherm(S(2));
isotherm{2} = isotherm2{1};
minlnP = [-Inf, minlnP2];
ads_pot{2} = [];
inv_ads_pot{2} = [];
[Q, x, err, lnP0, psi] = IAST_solve(M, S, 'mode', 2, 'ads_pot', ads_pot, 'inv_ads_pot', inv_ads_pot, 'options', options);

%% Plotting
semilogx(M(:,2),ethanol_water(:,2),'rs',M(:,2),ethanol_water(:,5),'bo',M(:,2),Q(:,2),'md',M(:,2),Q(:,1),'c^');
xlabel('$p$ [Pa]'); ylabel('$N$ [molec/uc]');
legend('sim\_EtOH','sim\_H2O','IAST\_EtOH','IAST\_H2O','Location','NorthEastOutside');

%% water/methanol/ethanol
S={[water(:,1),water(:,2)],[methanol(:,1),methanol(:,2)],[ethanol(:,1),ethanol(:,2)]};
M=[ternary(:,9)*1000,ternary(:,3)*1000,ternary(:,6)*1000];
[Q, x, err, lnP0, psi] = IAST_solve(M, S, 'mode', 1, 'options', options);

plot([0,1],[0,1],'k');hold on
plot(ternary(:,2)/(ternary(:,8)+ternary(:,2)+ternary(:,5)),Q(:,2)/(Q(:,1)+Q(:,2)+Q(:,3)),'mv');hold on
plot(ternary(:,5)/(ternary(:,8)+ternary(:,2)+ternary(:,5)),Q(:,3)/(Q(:,1)+Q(:,2)+Q(:,3)),'ch');
xlabel('$x_{\textrm{sim}}$'); ylabel('$x_{\textrm{IAST}}$');