%% IAST DRIVER: runs IAST calculations 
% see IAST_solve

%H2S/CO2/CH4
clc; clear all

function run_ternary(z)

    h2s=importdata(sprintf('data/ternary_%s/%s_H2S_tra_298K.csv', z, z)); % p/q pairs
    co2=importdata(sprintf('data/ternary_%s/%s_CO2_298K.csv', z, z));
    ch4=importdata(sprintf('data/ternary_%s/%s_methane_298K.csv', z, z));
    ternary_raw=importdata(sprintf('data/ternary_%s/%s_298K_3comp_mixture_226.csv', z, z));
    ternary=ternary_raw.data;

    options=optimset('Display','iter');
    S={[h2s.data(:,1),h2s.data(:,2)],[co2.data(:,1),co2.data(:,2)],[ch4.data(:,1),ch4.data(:,2)]};
    M=ternary(:,1)*[2,2,6]/10; %q_H2S, q_CO2, q_CH4
    [Q, x, err, lnP0, psi] = IAST_solve(M, S, 'mode', 1, 'options', options);

    run('startup.m')
    plot([0,1],[0,1],'k');hold on
    plot(ternary(:,3)/(ternary(:,2)+ternary(:,3)+ternary(:,4)),Q(:,2)/(Q(:,1)+Q(:,2)+Q(:,3)),'mv');hold on
    plot(ternary(:,4)/(ternary(:,2)+ternary(:,3)+ternary(:,4)),Q(:,3)/(Q(:,1)+Q(:,2)+Q(:,3)),'ch');
    xlabel('x_{sim}'); ylabel('x_{IAST}');

end


zeo='AEL-1';
t = timeit(@run_ternary(zeo))
