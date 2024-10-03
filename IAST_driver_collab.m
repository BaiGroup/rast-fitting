%% IAST DRIVER: runs IAST calculations 
% see IAST_solve

clc; clear all

function [Q, x, err, lnP0, psi] = run_ternary(M, S)
    
    options=optimset('Display','iter');
    M_pressure = [M{1}(:,1), M{2}(:,1), M{3}(:,1)]; %p_H2S, p_CO2, p_CH4
    [Q, x, err, lnP0, psi] = IAST_solve(M_pressure, S, 'mode', 1, 'options', options);

end

% choose zeolite + system
temp = 298;
zeo='AEL-1';
ratio = [2,2,6];
comp = {'H2S', 'CO2', 'methane'};

% load data
comp1=importdata(sprintf('data/ternary_%s/%s_%s_tra_%dK.csv', zeo, zeo, comp{1}, temp)); % p/q pairs
comp2=importdata(sprintf('data/ternary_%s/%s_%s_%dK.csv', zeo, zeo, comp{2}, temp));
comp3=importdata(sprintf('data/ternary_%s/%s_%s_%dK.csv', zeo, zeo, comp{3}, temp));
ternary_raw=importdata(sprintf('data/ternary_%s/%s_%dK_3comp_mixture_226.csv', zeo, zeo, temp));
ternary=ternary_raw.data;

S={[comp1.data(:,1), comp1.data(:,2)], [comp2.data(:,1), comp2.data(:,2)], ...
    [comp3.data(:,1),comp3.data(:,2)]}; %order: H2S, CO2, CH4
M={[ternary(:,1)*ratio(1)/10, ternary(:,2)], [ternary(:,1)*ratio(2)/10, ternary(:,3)], ...
    [ternary(:,1)*ratio(3)/10, ternary(:,4)]}; 

% % run IAST
[Q, x, err, lnP0, psi] = run_ternary(M, S);

% use timeit to get runtime (s)
% func = @() run_ternary(M, S);
% t = timeit(func, 1);
% fprintf('For zeolite %s, IAST takes %.2f s\n', zeo, t)

% plot
run('startup.m')
plot([0,1],[0,1],'k');hold on
plot(M{2}(:,2)/(M{1}(:,2)+M{2}(:,2)+M{3}(:,2)),Q(:,2)/(Q(:,1)+Q(:,2)+Q(:,3)),'mv');hold on
plot(M{3}(:,2)/(M{1}(:,2)+M{2}(:,2)+M{3}(:,2)),Q(:,3)/(Q(:,1)+Q(:,2)+Q(:,3)),'ch');
xlabel('x_{sim}'); ylabel('x_{IAST}');
