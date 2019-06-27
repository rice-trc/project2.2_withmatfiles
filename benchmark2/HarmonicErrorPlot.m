clear; clc;
close all; 
set(0,'defaultAxesTickLabelInterpreter', 'default');
set(0,'defaultTextInterpreter','latex'); 
set(0, 'DefaultLegendInterpreter', 'latex');

load('HarmonicErrorI.mat')

figure(1); set(gca, 'YScale', 'log')
hold on
semilogy(Harmonics,freq_err, 'LineWidth', 2)
semilogy(Harmonics,amp_err, 'LineWidth', 2)
semilogy(Harmonics,tolerence,'--k', 'LineWidth', 2)
dot = interp1(Harmonics,amp_err,9);
semilogy(9,dot,'ko', 'MarkerFaceColor', 'k')
dot = interp1(Harmonics,freq_err,9);
semilogy(9,dot,'ko', 'MarkerFaceColor', 'k')
hold off

errmax = 0.01*1e-2;
xlabel('Number of Harmonics')
ylabel('Relative Error (\%)')
legend('Resonant frequency', 'Resonant amplitude', sprintf('Threshold: %.2f \\%%', errmax*100), sprintf('Convergence: $N_h$=%d',9))

load('HarmonicErrorIII.mat')

figure(2); set(gca, 'YScale', 'log')
hold on
semilogy(Harmonics,freq_err, 'LineWidth', 2)
semilogy(Harmonics,amp_err, 'LineWidth', 2)
semilogy(Harmonics,tolerence,'--k', 'LineWidth', 2)
dot = interp1(Harmonics,amp_err,3);
semilogy(3,dot,'ko', 'MarkerFaceColor', 'k')
dot = interp1(Harmonics,freq_err,3);
semilogy(3,dot,'ko', 'MarkerFaceColor', 'k')
hold off

errmax = 0.01*1e-2;
xlabel('Number of Harmonics')
ylabel('Relative Error (\%)')
legend('Resonant frequency', 'Resonant amplitude', sprintf('Threshold: %.2f \\%%', errmax*100), sprintf('Convergence: $N_h$=%d',3))