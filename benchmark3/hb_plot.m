% close all
% clearvars

srcpath = '../src/nlvib/SRC';
addpath(srcpath);
srcpath = '../src/matlab';
addpath(genpath(srcpath));

savefig = false;

%% use Latex
set(0,'defaultAxesTickLabelInterpreter', 'default');
set(0,'defaultTextInterpreter','latex'); 
set(0, 'DefaultLegendInterpreter', 'latex'); 

% time points in period
Ntd = 1e3;

%% load data
% load data from mat-file if called from this script
% https://stackoverflow.com/a/35359201
stack=dbstack;
if numel(stack) == 1 || exist('frfdata','var') == 0
    frfdata = load('frf.mat');
    nmadata = load('nma.mat');
end

frf = process_frf(frfdata,Ntd);
nma = process_nma(nmadata,Ntd);

% % Illustrate frequency response
% figure; hold on; box on
% aa = gobjects(frf.nex,1);
% for j = 1:frf.nex
%     aa(j) = plot(frf.Om{j}/2/pi,frf.apeak{j}*1000,'k-','linewidth',1.5);
%     legend(aa(j), sprintf('F = %.2f N', frf.exc_lev(j)));
% end
% legend(aa, 'Location', 'northeast')
% xlabel('$f_{\mathrm{ex}}$ in Hz');
% ylabel('$\hat{w}_{L/2}$ in mm');
% title(['thickness = ' num2str(thickness*1000) 'mm']);



% % Determine total energy in the system from the displacement and velocity
% % at t=0
% energy = zeros(size(a_NMA));
% % loop over countinuation parameter(omega)
% for i=1:length(nma.Om)
%     Qi = reshape(Q_HB(:,i),n,2*H+1);
%     q0 = Qi(:,1) + sum(Qi(:,2:2:end),2);
%     u0 = sum(Qi(:,3:2:end),2)*om_HB(i);
%     energy(i) = 1/2*u0'*beam.M*u0 + 1/2*q0'*beam.K*q0 + knl*q0(n-1)^4/4;
% end

%% set default line style
% https://www.mathworks.com/help/matlab/graphics_transition/why-are-plot-lines-different-colors.html
set(groot,'defaultAxesColorOrder',[1 0 0;0 1 0;0 0 1;0 0.5 0;0.75 0 0.75],...
    'defaultAxesLineStyleOrder','--|:|-')

hfig = {};
%% FEP

% % Modal frequency vs. amplitude
% figure; hold on;
% plot(log10(energy),Om_nma/(2*pi),'k-o');
% xlabel('log10(energy)'); ylabel('modal frequency in Hz');
% % set(gca,'ylim',[20 50]);

%% NFRC

% peak amplitude
figure;
hold on
plot(nma.Om/2/pi,1000*abs(nma.apeak),'color',[1 .2 .3],'linewidth',2)

slabel = cell(1,frf.nex);
for j = 1:frf.nex
    plot(frf.Om{j}/2/pi,frf.apeak{j}*1000,'k-','linewidth',1.5);
    slabel{j} = sprintf('F = %0.2f N',frf.exc_lev(j));
end
legend(['nm - backbone', slabel], 'Location', 'northeast')

xlim([260 420])
xlabel('$f_{\mathrm{ex}}$ in Hz');
ylabel('$\hat{w}_{L/2}$ in mm');
title(['thickness = ' num2str(thickness*1000) 'mm']);
hfig(end+1) = {{gcf, "nfrc_peak"}};

% RMS amplitude
figure;
hold on
plot(nma.Om/2/pi,1000*abs(nma.arms),'color',[1 .2 .3],'linewidth',2)

slabel = cell(1,frf.nex);
for j = 1:frf.nex
    plot(frf.Om{j}/2/pi,frf.arms{j}*1000,'k-','linewidth',1.5);
    slabel{j} = sprintf('F = %0.2f N',frf.exc_lev(j));
end
legend(['nm - backbone', slabel], 'Location', 'northeast')

xlim([260 420])
xlabel('$f_{\mathrm{ex}}$ in Hz');
ylabel('$RMS \hat{w}_{L/2}$ in mm');
title(['thickness = ' num2str(thickness*1000) 'mm']);
hfig(end+1) = {{gcf, "nfrc_rms"}};

%% bending modes

figure;
plot(nma.Om/2/pi,nma.a_q,'linewidth',1.5)
legend('$q_1$ (first bending)','$q_2$ (third bending)',...
    '$q_3$ (fift bending)')
xlabel('$\omega$ in Hz');
ylabel('modal coordinates $q$ of NMA')
xlim([min(nma.Om/2/pi) max(nma.Om/2/pi)]);
title(['thickness = ' num2str(thickness*1000) 'mm']);
hfig(end+1) = {{gcf, "bending_modes"}};

%% save figures

if savefig
    path = 'fig/';
    for i=1:length(hfig)
        h = hfig{i}{1};
        fname = hfig{i}{2};
        % change figure background from standard grey to white
        set(h, 'Color', 'w');
        export_fig(h, strcat(path,fname), '-pdf', '-png');
    end
    
end

%% remove default linestyle
set(groot,'defaultAxesLineStyleOrder','remove')
set(groot,'defaultAxesColorOrder','remove')