clear; clc;
close all; 
set(0,'defaultAxesTickLabelInterpreter', 'default');
set(0,'defaultTextInterpreter','latex'); 
set(0, 'DefaultLegendInterpreter', 'latex'); 

load FlatBeam_FRF.mat

% Illustrate phase diagram
figure(1)
hold on
for iex=1:length(exc_lev)
    plot(Om{iex}/2/pi, phase{iex}, '.', 'LineWidth', 2);
end

for iex=1:length(exc_lev)
    res{iex} = interp1(phase{iex},Om{iex}/2/pi,-90);
    plot(res{iex},-90,'ok', 'LineWidth', 2)  % find resonance point
end

xlabel('Frequency (Hz)')
ylabel('Response Phase (degs)')
legend('Fex = 10','Fex = 40','Fex = 60','Fex = 80','Fex = 100',...
    'Location','E');
hold off

% Illustrate frequency response
figure (2); hold on; box on
for iex=1:length(exc_lev)
plot(Om{iex}/2/pi,a_w_L_2{iex}*1000,'linewidth',1.5);
end
plot(om_HB/2/pi,1000*abs(a_w_L_2_NMA),'--k','linewidth',2)
legend('Fex = 10','Fex = 40','Fex = 60','Fex = 80','Fex = 100','EPMC',...
    'Location','E');
for iex=1:length(exc_lev)
peak_amp = interp1(Om{iex}/2/pi,a_w_L_2{iex}*1000,res{iex}); % mark resonance point
plot(res{iex},peak_amp,'ok', 'LineWidth', 2)
end
xlabel('Omega')
ylabel('Amplitude (mm)')
title('Flat Beam FRF -- Mode 1')
xlim([100 650])
ylim([0 2])

plot(om_HB/2/pi,1000*abs(a_w_L_2_NMA),'--k','linewidth',2)
legend('Fex = 10','Fex = 40','Fex = 60','Fex = 80','Fex = 100','EPMC',...
    'Location','E');
hold off

figure (3)
Y_HB_1 = Q_HB(n+(1:n),:)-1i*Q_HB(2*n+(1:n),:);
% semilogx(abs(Y_HB_1(3,:)),del_HB*1e2,'b-','LineWidth',2);
semilogx(1000*abs(a_w_L_2_NMA),del_HB*1e2,'b-','LineWidth',2);
xlabel('amplitude in mm'); ylabel('modal damping ratio in %')
title('Damping Ratio')

figure (4)
semilogx(1000*abs(a_w_L_2_NMA),om_HB/2/pi,'linewidth',2)
xlabel('amplitude in mm'); ylabel('Frequency (Hz)')