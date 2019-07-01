clear; clc;
close all; 
set(0,'defaultAxesTickLabelInterpreter', 'default');
set(0,'defaultTextInterpreter','latex'); 
set(0, 'DefaultLegendInterpreter', 'latex'); 

%% FlatBeam_FRF_I.mat
load FlatBeam_FRF_I.mat

% Illustrate phase diagram
figure(1)
hold on
for iex=1:length(exc_lev)
    plot(Om{iex}/2/pi, phase{iex}, '-', 'LineWidth', 2);
end

for iex=1:length(exc_lev)
    res{iex} = interp1(phase{iex},Om{iex}/2/pi,-90);
    plot(res{iex},-90,'ok', 'LineWidth', 2)  % find resonance point
end

xlabel('Forcing Frequency $\omega$ (Hz)')
ylabel('Response phase (degs)')
legend('Fex = 10','Fex = 40','Fex = 60','Fex = 80','Fex = 100',...
    'Location','E');
hold off

% Illustrate frequency response
figure (2); hold on; box on
for iex=1:length(exc_lev)
plot(Om{iex}/2/pi,a_w_L_2{iex},'linewidth',1.5);
end
plot(om_HB/2/pi,abs(a_w_L_2_NMA),'--k','linewidth',2)
legend('Fex = 10','Fex = 40','Fex = 60','Fex = 80','Fex = 100','EPMC',...
    'Location','E');
for iex=1:length(exc_lev)
peak_amp = interp1(Om{iex}/2/pi,a_w_L_2{iex},res{iex}); % mark resonance point
plot(res{iex},peak_amp,'ok', 'LineWidth', 2)
end
xlabel('Forcing Frequency $\omega$ (Hz)')
ylabel('RMS Response Displacement Amplitude (m)')
% title('Flat Beam FRF -- Mode 1')
xlim([100 650])

plot(om_HB/2/pi,abs(a_w_L_2_NMA),'--k','linewidth',2)
legend('Fex = 10','Fex = 40','Fex = 60','Fex = 80','Fex = 100','EPMC',...
    'Location','E');
hold off

figure (3)
semilogx(a_NMA,del_HB*1e2,'b-','LineWidth',2);
% semilogx(1000*abs(a_w_L_2_NMA),del_HB*1e2,'b-','LineWidth',2);
xlabel('Modal Amplitude'); ylabel('Damping Factor')
title('Damping Ratio')

figure (4)
semilogx(a_NMA,om_HB/2/pi,'linewidth',2);
xlabel('Modal Amplitude'); ylabel('Natural Frequency (Hz)')

% figure (14)
% energy1 = sum(1/2*M*abs(Y_HB_1).^2.*om_HB.^2 + 1/2*K*abs(Y_HB_1).^2);
energy1 = a_NMA.^2.*om_HB.^2;
om1 = om_HB;
% semilogx(energy1,om_HB/2/pi,'linewidth',2);

%% FlatBeam_FRF_III.mat
load FlatBeam_FRF_III.mat

% Illustrate phase diagram
figure(5)
hold on
for iex=1:length(exc_lev)
    plot(Om{iex}/2/pi, phase{iex}, '-', 'LineWidth', 2);
end

for iex=1:length(exc_lev)
    res{iex} = interp1(phase{iex},Om{iex}/2/pi,90);
    plot(res{iex},90,'ok', 'LineWidth', 2)  % find resonance point
end

xlabel('Forcing Frequency $\omega$ (Hz)')
ylabel('Response phase (degs)')
legend('Fex = 10','Fex = 40','Fex = 60','Fex = 80','Fex = 100',...
    'Location','E');
hold off
xlim ([1300 1600])

% Illustrate frequency response
figure (6); hold on; box on
for iex=1:length(exc_lev)
plot(Om{iex}/2/pi,a_w_L_2{iex},'linewidth',1.5);
end
plot(om_HB/2/pi,abs(a_w_L_2_NMA),'--k','linewidth',2)
legend('Fex = 10','Fex = 40','Fex = 60','Fex = 80','Fex = 100','EPMC',...
    'Location','E');
for iex=1:length(exc_lev)
peak_amp = interp1(Om{iex}/2/pi,a_w_L_2{iex},res{iex}); % mark resonance point
plot(res{iex},peak_amp,'ok', 'LineWidth', 2)
end
xlabel('Forcing Frequency $\omega$ (Hz)')
ylabel('RMS Response Displacement Amplitude (m)')
% title('Flat Beam FRF -- Mode 1')
xlim ([1300 1600])
ylim ([0 2.5e-4])

plot(om_HB/2/pi,abs(a_w_L_2_NMA),'--k','linewidth',2)
legend('Fex = 10','Fex = 40','Fex = 60','Fex = 80','Fex = 100','EPMC',...
    'Location','E');
hold off

figure (7)
semilogx(a_NMA,del_HB*1e2,'b-','LineWidth',2);
% semilogx(1000*abs(a_w_L_2_NMA),del_HB*1e2,'b-','LineWidth',2);
xlabel('Modal Amplitude'); ylabel('Damping Factor')
title('Damping Ratio')

figure (8)
semilogx(a_NMA,om_HB/2/pi,'linewidth',2);
xlabel('Modal Amplitude'); ylabel('Natural Frequency (Hz)')

% figure (13)
% energy2 = a_NMA.^2.*om_HB.^2;
% om2 = om_HB;
% 
% energy = linspace(3e-4, 98, 10000);
% om11 = interp1(energy1,om1,energy);
% om22 = interp1(energy2,om2,energy);
% om_ratio = om22./om11;
% 
% semilogx(energy,om_ratio,'linewidth',2);

%% FlatBeam_FRF_V.mat
load FlatBeam_FRF_V.mat

% Illustrate phase diagram
figure(9)
hold on
for iex=1:length(exc_lev)
    plot(Om{iex}/2/pi, phase{iex}, '-', 'LineWidth', 2);
end

for iex=1:length(exc_lev)
    res{iex} = interp1(phase{iex},Om{iex}/2/pi,-90);
    plot(res{iex},-90,'ok', 'LineWidth', 2)  % find resonance point
end

xlabel('Forcing Frequency $\omega$ (Hz)')
ylabel('Response phase (degs)')
legend('Fex = 10','Fex = 40','Fex = 60','Fex = 80','Fex = 100',...
    'Location','E');
hold off
xlim ([3480 3600])

% Illustrate frequency response
figure (10); hold on; box on
for iex=1:length(exc_lev)
plot(Om{iex}/2/pi,a_w_L_2{iex},'linewidth',1.5);
end
plot(om_HB/2/pi,abs(a_w_L_2_NMA),'--k','linewidth',2)
legend('Fex = 10','Fex = 40','Fex = 60','Fex = 80','Fex = 100','EPMC',...
    'Location','E');
for iex=1:length(exc_lev)
peak_amp = interp1(Om{iex}/2/pi,a_w_L_2{iex},res{iex}); % mark resonance point
plot(res{iex},peak_amp,'ok', 'LineWidth', 2)
end
xlabel('Forcing Frequency $\omega$ (Hz)')
ylabel('RMS Response Displacement Amplitude (m)')
% title('Flat Beam FRF -- Mode 1')
xlim ([3480 3600])
ylim ([0 3.1e-5])

plot(om_HB/2/pi,abs(a_w_L_2_NMA),'--k','linewidth',2)
legend('Fex = 10','Fex = 40','Fex = 60','Fex = 80','Fex = 100','EPMC',...
    'Location','E');
hold off

figure (11)
Y_HB_1 = Q_HB(n+(1:n),:)-1i*Q_HB(2*n+(1:n),:);
semilogx(a_NMA,del_HB*1e2,'b-','LineWidth',2);
% semilogx(1000*abs(a_w_L_2_NMA),del_HB*1e2,'b-','LineWidth',2);
xlabel('Modal Amplitude'); ylabel('Damping Factor')
title('Damping Ratio')

figure (12)
semilogx(a_NMA,om_HB/2/pi,'linewidth',2);
xlabel('Modal Amplitude'); ylabel('Natural Frequency (Hz)')
