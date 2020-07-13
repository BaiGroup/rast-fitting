%% Load data
load('298K.mat');
set(0,'DefaultTextInterpreter','latex');

%% water/methanol
S={[water(:,1),water(:,2)],[methanol(:,1),methanol(:,2)]};
M=[methanol_water(:,6)*1000,methanol_water(:,3)*1000];
Q=IAST_solve(M,S);

semilogx(M(:,2),methanol_water(:,2),'rs',M(:,2),methanol_water(:,5),'bo',M(:,2),Q(:,2),'md',M(:,2),Q(:,1),'c^');
xlabel('$p$ [Pa]'); ylabel('$N$ [molec/uc]');
legend('sim\_MeOH','sim\_H2O','IAST\_MeOH','IAST\_H2O','Location','NorthEastOutside');

%% water/ethanol
S={[water(:,1),water(:,2)],[ethanol(:,1),ethanol(:,2)]};
M=[ethanol_water(:,6)*1000,ethanol_water(:,3)*1000];
Q=IAST_solve(M,S);

semilogx(M(:,2),ethanol_water(:,2),'rs',M(:,2),ethanol_water(:,5),'bo',M(:,2),Q(:,2),'md',M(:,2),Q(:,1),'c^');
xlabel('$p$ [Pa]'); ylabel('$N$ [molec/uc]');
legend('sim\_EtOH','sim\_H2O','IAST\_EtOH','IAST\_H2O','Location','NorthEastOutside');

%% water/methanol/ethanol
S={[water(:,1),water(:,2)],[methanol(:,1),methanol(:,2)],[ethanol(:,1),ethanol(:,2)]};
M=[ternary(:,9)*1000,ternary(:,3)*1000,ternary(:,6)*1000];
Q=IAST_solve(M,S);

plot([0,1],[0,1],'k');hold on
plot(ternary(:,2)/(ternary(:,8)+ternary(:,2)+ternary(:,5)),Q(:,2)/(Q(:,1)+Q(:,2)+Q(:,3)),'mv');hold on
plot(ternary(:,5)/(ternary(:,8)+ternary(:,2)+ternary(:,5)),Q(:,3)/(Q(:,1)+Q(:,2)+Q(:,3)),'ch');
xlabel('$x_{\textrm{sim}}$'); ylabel('$x_{\textrm{IAST}}$');