% close all
% clearvars

srcpath = '../src/nlvib/SRC';
addpath(srcpath);
srcpath = '../src/matlab';
addpath(genpath(srcpath));

savefig = true;
hfig = {};
imod = [1];


%% use Latex
set(0,'defaultAxesTickLabelInterpreter', 'default');
set(0,'defaultTextInterpreter','latex'); 
set(0, 'DefaultLegendInterpreter', 'latex'); 

%% FRF, RMS amplitude

xlms = [260 420;];
Ntd = 1024;
for i=imod

frfdata = load(sprintf('data/frf_M%d.mat',i));
nmadata = load(sprintf('data/nma_M%d.mat',i));

frf = process_frf(frfdata,Ntd);
nma = process_nma(nmadata,Ntd);
cmap = distinguishable_colors(frf.nex,'k');

figure;
hold on
plot(nma.Om/2/pi,1000*abs(nma.arms),'k-','linewidth',1.5)

slabel = cell(1,frf.nex);
for j = 1:frf.nex
    plot(frf.Om{j}/2/pi,frf.arms{j}*1000,'-','Color',cmap(j,:),'linewidth',1.5);
    slabel{j} = sprintf('F = %0.2f N',frf.exc_lev(j));
end
legend(['nm - backbone', slabel], 'Location', 'nw')

xlim(xlms(i,:))
xlabel('$f_{\mathrm{ex}}$ in Hz');
ylabel('$RMS \hat{w}_{L/2}$ in mm');
title(['thickness = ' num2str(frf.thickness*1000) 'mm']);
hfig(end+1) = {{gcf, sprintf("nfrc_rms_M%d",i)}};

%% phase plot
figure;
hold on

slabel = cell(1,frf.nex);
for j = 1:frf.nex
    plot(frf.Om{j}/2/pi,frf.phase{j}(end,:),'-','Color',cmap(j,:),'linewidth',1.5);
    slabel{j} = sprintf('F = %0.2f N',frf.exc_lev(j));
end
legend(slabel, 'Location', 'nw')

xlim(xlms(i,:))
xlabel('$f_{\mathrm{ex}}$ in Hz');
ylabel('Response Phase (degs)')
title(['thickness = ' num2str(frf.thickness*1000) 'mm']);
hfig(end+1) = {{gcf, sprintf("nfrc_phase_M%d",i)}};

%% bending modes

figure;
plot(nma.Om/2/pi,nma.Qrms,'linewidth',1.5)
legend('$q_1$ (first bending)','$q_2$ (third bending)',...
    '$q_3$ (fift bending)')
xlabel('$\omega$ in Hz');
ylabel('modal coordinates $q$ of NMA')
xlim([min(nma.Om/2/pi) max(nma.Om/2/pi)]);
title(['thickness = ' num2str(frf.thickness*1000) 'mm']);
hfig(end+1) = {{gcf, "bending_modes"}};

end

%% NM freq-amp / freq-damp

cmap = distinguishable_colors(length(imod),'k');

figure;
hfig(end+1) = {{gcf, "nm_freq"}};
h1=axes;
hold on
xlabel('Modal amplitude (m)');
ylabel('Frequency (Hz)');

figure;
hfig(end+1) = {{gcf, "nm_amp"}};
h2=axes;
hold on
xlabel('Modal amplitude (m)');
ylabel('damping (\%)');

a1 = gobjects(length(imod));
a2 = gobjects(length(imod));
for i=imod
nmadata = load(sprintf('data/nma_M%d.mat',i));
nma = process_nma(nmadata,Ntd);

a1(i) = plot(h1, abs(nma.arms),nma.Om/2/pi,'Color', cmap(i,:),'linewidth',1.5);
a2(i) = plot(h2, abs(nma.arms),nma.xi*1e2,'Color', cmap(i,:),'linewidth',1.5);
legend(a1(i), sprintf('mode %d',i))
legend(a2(i), sprintf('mode %d',i))
end

set(h1, 'XScale', 'log')
set(h2, 'XScale', 'log')

%% HB convergence

cmap = distinguishable_colors(length(imod));

figure;
hold on
for i=imod
    hbc = load(sprintf('data/HConvDat_M%d.mat',i));
    
%     popt = {'Color',cmap(imod,:),'LineWidth', 2};
    h(i) = plot(hbc.Hs(1:end-2), hbc.Errs(1:end-2,1)*100, '-', 'Color',cmap(i,:),'LineWidth', 2);
    plot(hbc.Hs(1:end-2), hbc.Errs(1:end-2,2)*100, ':', 'Color',cmap(i,:),'LineWidth', 2)

    plot(hbc.Hs(hbc.ihconv), hbc.Errs(hbc.ihconv,:)*100, 'o',...
        'Color',cmap(i,:),'MarkerFaceColor', cmap(i,:))
    hline(hbc.errmax*100,'k--');  % treshold
    str{i} = sprintf('mode = %d', i);
end

% create lines for legend
L(1) = plot(nan, nan, 'k-','LineWidth', 2);
L(2) = plot(nan, nan, 'k:','LineWidth', 2);
L(3) = plot(nan, nan, 'ko','MarkerFaceColor','k');
L(4) = plot(nan, nan, 'k--');

set(gca, 'YScale', 'log')
xlabel('Number of Harmonics')
ylabel('Relative Error (\%)')
legend([h,L],{str{:},'Resonant frequency', 'Resonant amplitude',...
    'Convergence',sprintf('Threshold: %.2f \\%%', hbc.errmax*100) } )

%     ... 
%     sprintf('Convergence: $N_h$=%d', hbc.Nhconv),...
%     sprintf('Threshold: %.2f \\%%', hbc.errmax*100)} )

hfig(end+1) = {{gcf, "convergence"}};


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