clear; clc;
close all; 
set(0,'defaultAxesTickLabelInterpreter', 'default');
set(0,'defaultTextInterpreter','latex'); 
set(0, 'DefaultLegendInterpreter', 'latex');

load('HarmonicError.mat')

figure(1); set(gca, 'YScale', 'log')
hold on
semilogy(Harmonics,freq_err, 'LineWidth', 2)
semilogy(Harmonics,amp_err, 'LineWidth', 2)
semilogy(Harmonics,tolerence,'--k', 'LineWidth', 2)
hold off

xlabel('Number of Harmonics')
ylabel('Relative Error (%)')
legend('Resonance Frequency','Resonance Amplitude','Tolerence')
savefig('HarmonicError.fig')