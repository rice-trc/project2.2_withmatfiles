% estimate PNLSS model
%
%  1) estimate a nonparametric linear model from the input and 
%     noisy output data
%  2) estimate a parametric linear state-space model on the
%     nonparametric model
%  3) estimate the parameters in the full PNLSS model

close all
clearvars
% clc

srcpath = '../src/pnlss';
addpath(genpath(srcpath));
figpath = './fig/';

dataname = 'ms_full';
A = 15;

%%
% addnoise = true;
addnoise = false;
savefig = true;

load(sprintf('data/A%d_%s',A,dataname))

[Nt,P,R,p] = size(y);
[Nt,P,R,m] = size(u);

% convert to measurement points using modal shapes from Phi.
nPHIS = size(PHIS,1);
yres = zeros(Nt,P,R,nPHIS); 
for i = 1:nPHIS
    yres(:,:,:,i) = reshape(PHIS(i,:)*reshape(y,[],p)', [Nt,P,R,1]);
end
y = yres;
p = nPHIS;

% % use middle deflection
% y = PHI_L2*reshape(y,[],n)';
% y = reshape(y,[Nt,P,R]);


%% Add colored noise to the output
rng(10);
if addnoise
    noise = 1e-3*std(y(:,end,end))*randn(size(y)); % Output noise signal
    % Do some filtering
    noise(1:end-1,:,:) = noise(1:end-1,:,:) + noise(2:end,:,:);
    y = y + noise;
end

%% Separate the data in estimation, validation, and test set
% Last realization, last period for performance testing
utest = u(:,end,R,:); utest = reshape(utest,[],m);
ytest = y(:,end,R,:); ytest = reshape(ytest,[],p);

% One but last realization, last period for validation and model selection
uval = u(:,end,R-1,:); uval = reshape(uval,[],m);
yval = y(:,end,R-1,:); yval = reshape(yval,[],p);

% All other realizations for estimation. But remember to remove transient!
R = R-2;
Ptr = 3;
P = P-Ptr;
uest = u(:,Ptr:end,1:R,:);
yest = y(:,Ptr:end,1:R,:);

% standard deviation of the generated signal
uStd = mean(mean(std(uest)));

% lines = cell2mat(cellfun(@(c) c.lines, MS, 'UniformOutput', false));
lines = MS{1}.lines;
lines(lines==1) = [];
%% Estimate nonparametric linear model (BLA)

% m: number of inputs, p: number of outputs
uest = permute(uest,[1,4,3,2]); % N x m x R x P
yest = permute(yest,[1,4,3,2]); % N x p x R x P
covY = fCovarY(yest); % Noise covariance (frequency domain)


U = fft(uest); U = U(lines,:,:,:); % Input spectrum at excited lines
Y = fft(yest); Y = Y(lines,:,:,:); % Output spectrum at excited lines

% Estimate best linear approximation, total distortion, and noise distortion
% total and noise distortion averaged over P periods and R realizations
% total distortion level includes nonlinear and noise distortion
% G: FRF; covGML: noise + NL; covGn: noise (all only on excited lines)
[G,covGML,covGn] = fCovarFrf(U,Y); 
% figure;
% subplot(2,1,1); semilogy(freq(lines),abs(squeeze(G(:))));
% subplot(2,1,2); plot(freq(lines),rad2deg(angle(G(:))))
%% Estimate linear state-space model (frequency domain subspace)

% Choose model order
na   = 4;%[2,3,4];
maxr = 5;%20;
% Excited frequencies (normalized)
freq_norm = (lines-1)/Nt;

% Uncomment for uniform weighting (= no weighting)
% covGML = repmat(eye(1),[1 1 length(lines)]);

models = fLoopSubSpace(freq_norm,G,covGML,na,maxr,100);

% Extract linear state-space matrices from best model on validation data
Nval = length(uval);
tval = (0:Nval-1)/fs;
min_err = Inf;
min_na = NaN;
for n = na
    model = models{n};
    A = model{1}; B = model{2}; C = model{3}; D = model{4};
    [A,B,C] = dbalreal(A,B,C); % Compute balanced realization
    yval_hat = lsim(ss(A,B,C,D,1/fs),uval,tval);
    err = yval - yval_hat; 
    % Rms value of the last period of the error signal
    err = sqrt(mean(err(end-Nt+1:end).^2));
    if err < min_err
        min_na = n;
        min_err = err;
    end
end
% select the best model
model = models{min_na};
[A,B,C,D] = model{:};
% Balanced realization
[A,B,C] = dbalreal(A,B,C);


%% Estimate PNLSS model

% Average over periods (be careful that the data are truly steady state)
uest = mean(uest,4); 
yest = mean(yest,4);  % N x p x R
uest = permute(uest,[1,3,2]); % N x R x m
yest = permute(yest,[1,3,2]); % N x R x p
uest = reshape(uest,Nt*R,m); % Concatenate the data: N*P*R x m
yest = reshape(yest,Nt*R,p); % N*P*R x p

% Transient settings
NTrans = 2*Nt; % Add one period before the start of each realization
% Number of transient samples and starting indices of each realization
T1 = [NTrans 1+(0:Nt:(R-1)*Nt)]; 
T2 = 0; % No non-periodic transient handling

% Nonlinear terms
nx = [3];
ny = [];
whichtermsx = 'statesonly';
% whichtermsx = 'full';
whichtermsy = 'empty';

% Settings Levenberg-Marquardt optimization
MaxCount = 100;
lambda = 100;

% Choose model order
n = min_na;

% Initial linear model in PNLSS form
model = fCreateNLSSmodel(A,B,C,D,nx,ny,T1,T2);

% Set which monomials will be optimized
model.xactive = fSelectActive(whichtermsx,n,m,n,nx);
model.yactive = fSelectActive(whichtermsy,n,m,p,ny);

% nonlinear powers
% tmp=kron([1;1],model.xpowers); tmp(model.xactive,:)

% Output of the initial linear model on the estimation data
modellinest = model;
y_lin = fFilterNLSS(modellinest,uest);
errest_lin = yest-y_lin;

% Compute output error on validation data
% Change transient parameter for linear model on validation data
modellinval = model;modellinval.T1 = Nt;
yval_lin = fFilterNLSS(modellinval,uval);
errval_lin = yval-yval_lin;

%Compute output error on test data
modellintest = model;modellintest.T1 = 0;
ytest_lin = fFilterNLSS(modellintest,utest);
errtest_lin = ytest-ytest_lin;

% We do not use weighting
% for kk = size(covY,3):-1:1
%     W(:,:,kk) = fSqrtInverse(covY(:,:,kk)); % Frequency weighting
% end
W = [];

% Levenberg-Marquardt optimization
[~, y_mod, models_pnlss] = fLMnlssWeighted(uest,yest,model,MaxCount,W,lambda);
errest_nl = yest-y_mod;

% Search best model over the optimisation path on a fresh set of data
valerrs = zeros(length(models_pnlss),1);
valerr = cell(length(models_pnlss),1);
for i = 1:length(models_pnlss)
    models_pnlss(i).T1 = NTrans;
    yval_mod = fFilterNLSS(models_pnlss(i),uval); 
    valerr{i} = rms(yval - yval_mod,1);
    valerrs(i) = sum(valerr{i});
end

% Select the best model on the validation data to avoid overfitting
[min_valerr,i] = min(valerrs);
model = models_pnlss(i);

% Compute output error on validation data
modelval = model;modelval.T1 = NTrans;
yval_nl = fFilterNLSS(modelval,uval);
errval_nl = yval-yval_nl;
% Compute output error on test data
modeltest = model;modeltest.T1 = 0;
ytest_nl = fFilterNLSS(modeltest,utest);
errtest_nl = ytest-ytest_nl;

% error
err = [rms(errest_lin), rms(errest_nl); rms(errval_lin), rms(errval_nl);...
    rms(errtest_lin), rms(errtest_nl)];

fprintf('############# RMS errors #############\n')
fprintf('e_est_lin:\t %0.3e\t e_est_nl:\t %0.3e\n', err(1,:))
fprintf('e_val_lin:\t %0.3e\t e_val_nl:\t %0.3e\n', err(2,:))
fprintf('e_test_lin:\t %0.3e\t e_test_nl:\t %0.3e\n',err(3,:))


save('./data/pnlssout_try0.mat', 'modellinest', 'model','valerr','valerrs')

%% Results

% convert to continious time for modal analysis of linear part
sys = d2c(ss(model.A,model.B,model.C,model.D,1/fs));
sd = modal(sys.A,sys,C);
disp(['Nat freq' sd.wn newline 'damping' sd.zeta])

% similarity transform to get state space matrices in physical coordinates
[Ap,Bp,Cp,T] = ss2phys(sys.A,sys.B,sys.C);
sys_phys = ss(Ap,Bp,Cp,model.D);

% handle to figs. For saving plot
hfig = {};

%% optimization path / convergence
figure;
plot(db(valerrs));
hold on
plot(i,db(min_valerr),'r.','Markersize',10)
xlabel('Successful iteration number')
ylabel('Validation error [dB]')
title('Selection of the best model on a separate data set')
hfig(end+1) = {{gcf, 'convergence'}};

%% Estimation data
cmap = distinguishable_colors(3,'k');

plottime = [yest errest_lin errest_nl];
plottime = reshape(plottime, [Nt*R,p,3]);
plotfreq = fft(reshape(plottime,[Nt,R,p,3]));
plotfreq = squeeze(mean(plotfreq,2));

% plot time series individual
for i=1:p
figure;
plot(t(1:Nt*R),squeeze(plottime(:,i,:)))
xlabel('Time (s)')
ylabel('Output (errors)')
legend('output','linear error','PNLSS error')
errstr = sprintf('rms(y-y_mod) = %0.3e\n(= Output error of the best PNLSS model on the estimation data)\n',rms(yest(:,i)-y_mod(:,i)) );
title(sprintf('Estimation results p: %d\n%s',errstr))
disp(' ')
disp(errstr)
% disp(['rms(noise(:))/sqrt(P) = ' num2str(rms(noise(:))/sqrt(P)) ' (= Noise level)'])
hfig(end+1) = {{gcf, sprintf('est_time_p%d',p)}};
end

% plot time series together
figure;
plot(t(1:Nt*R),reshape(plottime, [Nt*R, p*3]))
xlabel('Time (s)')
ylabel('Output (errors)')
legend('output','linear error','PNLSS error')
errstr = sprintf('sum(rms(y-y_mod)) = %0.3e\n(= Output error of the best PNLSS model on the estimation data)\n',sum(rms(yest-y_mod,1)) );
title(['Estimation results' newline errstr])
disp(' ')
disp(errstr)
hfig(end+1) = {{gcf, sprintf('est_time')}};

% error in frequency domain
for i=1:p
figure;
hold on
plot(freq(1:end/2),db(squeeze(plotfreq(1:end/2,i,:))),'.') %,'Color',cmap)
plot(freq(1:end/2),db(sqrt(P*squeeze(covY(i,i,:)))), 'k.')
xlabel('Frequency (Hz)')
ylabel('Output (errors) (dB)')
legend('Output','Linear error','PNLSS error','noise')
title(sprintf('Estimation results p: %d',i))
hfig(end+1) = {{gcf, sprintf('est_err_p%d',i)}};
end

%% Validation data
plottime = [yval errval_lin errval_nl];
plottime = reshape(plottime, [Nt,p,3]);
plotfreq = fft(plottime);

for i=1:p
figure;
plot(freq(1:end/2),db(squeeze(plotfreq(1:end/2,i,:))),'.') % ,'Color',cmap(i,:))
xlabel('Frequency (Hz)')
ylabel('Output (errors) (dB)')
legend('Output','Linear error','PNLSS error')
title(sprintf('Validation results p: %d',i))
hfig(end+1) = {{gcf, sprintf('val_err_p%d',i)}};
end

%% Test, ie. newer seen data
plottime = [ytest errtest_lin errtest_nl];
plottime = reshape(plottime, [Nt,p,3]);
plotfreq = fft(plottime);

for i=1:p
figure;
plot(freq(1:end/2),db(squeeze(plotfreq(1:end/2,i,:))),'.')
xlabel('Frequency');
ylabel('Output (errors) (dB)');
legend('Output','Linear error','PNLSS error')
title(sprintf('Test results p: %d',i))
hfig(end+1) = {{gcf, sprintf('test_err_p%d',i)}};
end

%% BLA plot. We can estimate nonlinear distortion
% total and noise distortion averaged over P periods and M realizations
% total distortion level includes nonlinear and noise distortion

for i=1:p
figure;
hold on;
plot(freq(lines),db(abs(squeeze(G(i,1,:)))))
plot(freq(lines),db(abs(squeeze(covGn(i,i,:)))*R,'power'),'s')
plot(freq(lines),db(abs(squeeze(covGML(i,i,:)))*R,'power'),'*')
xlabel('frequency (Hz)')
ylabel('magnitude (dB)')
title(sprintf('Estimated BLA p: %d)',i))
legend('BLA FRF','Noise Distortion','Total Distortion','Location','nw')
hfig(end+1) = {{gcf, sprintf('bla_p%d',i)}};
end
%% save figures

if savefig
    path = './fig/';
    for i=1:length(hfig)
        h = hfig{i}{1};
        fname = hfig{i}{2};
        % change figure background from standard grey to white
        set(h, 'Color', 'w');
        export_fig(h, strcat(path,fname), '-pdf', '-png');
    end
    
end

