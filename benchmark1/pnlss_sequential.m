clc
clear all

srcpath = '../src/pnlss';
addpath(genpath(srcpath));
addpath('../src/matlab/')
figpath = './fig';

%% 
addnoise = false;
savefig  = true;

Alevels = [10 25 35 50 100 500 1000 10000 100000];
fs = 4096;

%% BLA with lowest forcing
load(sprintf('./data/ode45_multisine_A%d_F%d.mat',Alevels(1), fs))

freq = (0:Nt-1)*f0;
[Nt, P, R, n] = size(y);

% Use middle deflection
y = PHI_L2*reshape(y, [], n)';
y = reshape(y, [Nt, P, R]);

% %% Add colored noise to the output if required
rng(10)
if addnoise
    noise = 1e-3 *std(y(:,end,end))*randn(size(y));
    
    noise(1:end-1,:,:) = noise(1:end-1,:,:) + noise(2:end,:,:);
    y = y+noise;
end

% Separate the data in estimation, validation and test set
utest = u(:, end, end);  utest = utest(:);
ytest = y(:, end, end);  ytest = ytest(:);

uval = u(:, end, R-1);  uval = uval(:);
yval = y(:, end, R-1);  yval = yval(:);

% All other repeats for estimation
R = R-2;
Ptr = 6;
P = P-Ptr;
u = u(:, Ptr:end, 1:R);
y = y(:, Ptr:end, 1:R);

% Standard deviation of generated signal
uStd = mean(mean(std(u)));

% Non-parametric BLA model

% m: number of inputs, p: number of outputs
u = permute(u, [1,4,3,2]); % N x m x R x P
y = permute(y, [1,4,3,2]); % N x p x R x P
covY = fCovarY(y);  % Noise covariance (frequency domain)

lines = (f1/f0+1):(f2/f0+1);
U = fft(u);  U = U(lines, :, :, :);  % Input Spectrum at excited lines
Y = fft(y);  Y = Y(lines, :, :, :);  % Output Spectrum at excited lines

[G, covGML, covGn] = fCovarFrf(U, Y);

% Estimate linear state-space model (frequency domain subspace)
% Model order
na = 2;
maxr = 20;
% Excited frequencies (normed)
freqs_norm = (lines-1)/Nt;

models = fLoopSubSpace(freqs_norm, G, covGML, na, maxr, 100);

% Extract linear state-space matrices from best model on validation data
Nval = length(uval);
tval = (0:Nval-1)/fs;
min_err = Inf;
min_na = NaN;
for n = na
    model = models{n};
    A = model{1}; B = model{2}; C = model{3}; D = model{4};
    [A,B,C] = dbalreal(A,B,C);  % Balance realizations
    yval_hat = lsim(ss(A,B,C,D,1/fs), uval, tval);

    err = yval - yval_hat; 
    err = sqrt(mean(err(end-Nt+1:end).^2));
    if err < min_err
        min_na = n;
        min_err = err;
    end
end
min_na = na(1);
% Select the best model
model = models{min_na};
[A,B,C,D] = model{:};
[A,B,C] = dbalreal(A,B,C);

%% Misc Settings
% Transients
Ntrans = Nt;
T1 = [Ntrans 1+(0:Nt:(R-1)*Nt)];
T2 = 0;

m = size(u,2);
p = size(y,2);

nx = [3];
ny = [];
whichtermsx = 'statesonly';
whichtermsy = 'empty';

% Optimization
MaxCount = 100;
lambda = 100;

% model order 
n = min_na;

% BLA models with zero nonlinear coefficients
model = fCreateNLSSmodel(A, B, C, D, nx, ny, T1, T2);
model.xactive = fSelectActive(whichtermsx, n, m, n, nx);
model.yactive = fSelectActive(whichtermsy, n, m, p, ny);

modellinest = model;
modellinval = model;  modellinval.T1 = Nt;
modellintest = model; modellintest.T1 = 0;

% Weighting
W = [];  

%% Sequential PNLSS
errormeasures = cell(size(Alevels));
seqmodels = cell(size(Alevels));
modelguess = modellintest;
for ia=1:length(Alevels)
    load(sprintf('./data/ode45_multisine_A%d_F%d.mat',Alevels(ia), fs), 'u', 'y');
    
    % Kick off transients
    u = u(:, Ptr:end, 1:R);
    y = y(:, Ptr:end, 1:R);
    
    % Separate data
    utest = u(:, end, end);  utest = utest(:);
    ytest = y(:, end, end);  ytest = ytest(:);

    uval = u(:, end, R-1);  uval = uval(:);
    yval = y(:, end, R-1);  yval = yval(:);
    
    u = permute(u, [1,4,3,2]);  % N x m x R x P
    y = permute(y, [1,4,3,2]);  % N x p x R x P
    covY = fCovarY(y);  % Noise covariance (Frequency domain)
    
    u = mean(u, 4);
    y = mean(y, 4);
    m = size(u, 2);
    p = size(y, 2);
    uc = u(:);
    yc = y(:);
    
    % Simulate others
    y_linest = fFilterNLSS(modellinest, uc);
    err_linest = yc-y_linest;
    
    y_linval = fFilterNLSS(modellinval, uc);
    err_linval = yc-y_linval;
    
    y_lintest = fFilterNLSS(modellintest, uc);
    err_lintest = yc-y_lintest;
    
    y_nlinit = fFilterNLSS(model, uc);
    err_nlinit = yc-y_nlinit;

    % PNLSS optimization
    [~, y_mod, models_pnlss] = fLMnlssWeighted(uc, yc, modelguess, MaxCount, W, lambda);
 	[~, y_mod, models_pnlss] = fLMnlssWeighted(uc, yc, models_pnlss(end-5), MaxCount, W, lambda);
    err_nlest = yc-y_mod;
    
    % Choose best model from PNLSS (using validation data)
    valerrs = [];
    for i=1:length(models_pnlss)
        models_pnlss(i).T1 = Ntrans;
        yval_mod = fFilterNLSS(models_pnlss(i), uval);
        valerr = yval - yval_mod;
        valerrs = [valerrs; rms(valerr)];
    end
    
    [min_valerr, i] = min(valerrs);
    model = models_pnlss(i);    % Updated PNLSS model
    modelguess.A = model.A;
    modelguess.B = model.B;
    modelguess.C = model.C;
    modelguess.D = model.D;
    modelguess.E = model.E;
    modelguess.F = model.F;

    % Compute output error on validation data
    modelval = model;  modelval.T1 = Ntrans;
    y_nlval = fFilterNLSS(modelval, uval);
    err_nlval = yval-y_nlval;
    % Compute output error on test data
    modeltest = model;  modeltest.T1 = Ntrans;
    y_nltest = fFilterNLSS(modeltest, utest);
    err_nltest = ytest-y_nltest;
    
    err = ...
        [rms(err_linest), rms(err_nlest); ...
        rms(err_linval), rms(err_nlval);...
        rms(err_lintest), rms(err_nltest)];
    
    errormeasures{ia} = err;
    seqmodels{ia} = model;
    
    save(sprintf('./data/pnlssmodel_A%d_F%d.mat', Alevels(ia), fs), 'model', 'err')
    
    fprintf('Done %d/%d\n', ia, length(Alevels))
end

save('./data/seqpnlssmodels.mat', 'Alevels', 'errormeasures', 'seqmodels');