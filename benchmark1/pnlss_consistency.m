clc
clear all

%%
% famp = 150;
% linmodel   = load(sprintf('./data/ode45_multisine_f%d.mat',famp), 'model', 'fs', 'A');
% pnlssmodel = load(sprintf('./data/pnlssout_f%d.mat',famp), 'model');

linmodel   = load('./data/ode45_multisine.mat', 'model', 'fs');
pnlssmodel = load('./data/pnlssout_try0.mat', 'model');

lm = ss(linmodel.model.A, linmodel.model.B, linmodel.model.C, linmodel.model.D)
pm = d2c(ss(pnlssmodel.model.A, pnlssmodel.model.B, pnlssmodel.model.C, pnlssmodel.model.D, 1/linmodel.fs));

[Ap, Bp, Cp, Tp] = ss2phys(pm.A, pm.B, pm.C);
% Bp = Bp/linmodel.A;
% Cp = Cp*linmodel.A;

pmphy = ss(Ap, Bp, Cp, pm.D)
Bp(1)
Bp(2)

%%
Xpowers = pnlssmodel.model.xpowers(:,1:2);

linmodel.E = zeros(2, size(Xpowers,1));
linmodel.E(2, find(Xpowers(:,1)==3 & Xpowers(:,2)==0)) = linmodel.model.nlcof.coef;
% linmodel.Tmat = NLCOEF_TFMMATS(linmodel.Phi, Xpowers);

Tfmat = inv(Tp);
pnlssmodel.Tmat = NLCOEF_TFMMATS(Tp, Xpowers);

Ep = inv(Tfmat)*pnlssmodel.model.E*pnlssmodel.Tmat*linmodel.fs
El = linmodel.E
