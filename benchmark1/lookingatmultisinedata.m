clc
clear all

%% Histograms
Alevels=[10 50 100 1000 10000 100000];
fs = 4096;
sp = [3 2];

figure(1)
clf()

figure(2)
clf()
CFs = zeros(size(Alevels));
SPs = zeros(size(Alevels));

fprintf('======================\n');
fprintf('Signal Quality Metrics\n');
fprintf('----------------------\n');
fprintf('Alevel     CF  SD/Peak\n');
fprintf('----------------------\n');
for ia=1:length(Alevels)
    load(sprintf('./data/ode45_multisine_A%d_F%d.mat',Alevels(ia),fs), 'u', 'y');
    
    figure(1)
    subplot(sp(1), sp(2), ia)
    histogram(u(:,1,1), 'Normalization', 'pdf');
    xlabel('Force amplitude (N)')
    ylabel('pdf')
    title(sprintf('A = %d N', Alevels(ia)))    

    figure(2)
    subplot(sp(1), sp(2), ia)
    histogram(y(:,1,1), 'Normalization', 'pdf'); 
    xlabel('Displacement')
    ylabel('pdf')
    title(sprintf('A = %d N', Alevels(ia)))

    CFs(ia) = max(u(:,1,1))/rms(u(:,1,1));
    SPs(ia) = sqrt(var(u(:,1,1)))/max(u(:,1,1));
    
    fprintf('%d     %f  %f\n', Alevels(ia), CFs(ia), SPs(ia));
end
fprintf('======================\n');

figure(1)
print('./fig/Excitation_Hists.eps', '-depsc')
pause(1)

figure(2)
print('./fig/Response_Hists.eps', '-depsc')
%% Phase-Space

figure(3)
clf()
Alevels=[10 50 100 1000 10000 100000];
fs = 4096;
sp = [3 2];
for i=1:length(Alevels)
    load(sprintf('./data/ode45_multisine_A%d_F%d.mat',Alevels(i),fs), 'y', 'ydot');
    
    subplot(sp(1), sp(2), i)
%     plot(y(:,1,1), ydot(:,1,1), '.', 'MarkerSize', 0.1)
    scatter_kde(y(:,1,1), ydot(:,1,1), '.', 'MarkerSize', 0.1)
    yy=colorbar();
    ylabel(yy, 'pde')
    xlabel('y')
    ylabel('dy/dt')
    title(sprintf('A = %d N', Alevels(i)))
end
print('./fig/MS_stateplanekde.eps', '-depsc')