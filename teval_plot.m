%-------------------------------------------------------------------------%
% teval_plot - Jeremy Turner
% 
% teval_plot.m will plot the performance data of the four methods, assuming
% that each method has been analyzed by teval.m
% ------------------------------------------------------------------------%

clear all

f = {'DAE Linkage_teval.mat', 'DAE Pendulum_teval.mat', 'Lagrange Pendulum_teval.mat', 'Newton-Euler Pendulum_teval.mat'};

nm = {};
tavm = {};
tsdm = {};
ms = zeros(1, 4);
fm = {};

figure(1)
hold on
figure(2)
hold on
for i=1:length(f)
    
    load(f{i});
    
    nm{i} = ns;
    tavm{i} = tav;
    tsdm{i} = tsd;
    ms(i) = m;
    fm{i} = fs;
    
    figure(1)
    errorbar(ns, tav, tsd); xlabel('n'); ylabel('t [s]');
    figure(2)
    errorbar(find(ns<=10), tav(find(ns<=10)), tsd(find(ns<=10))); xlabel('n'); ylabel('t [s]');
    
end
figure(1)
legend(fm);
title('Method performance vs. Number of Links')
figure(2)
legend(fm);
title('Method performance vs. Number of Links')