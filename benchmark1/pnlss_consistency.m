clc
clear all

%%
linmodel   = load('./data/ode45_multisine.mat', 'model', 'fs');
pnlssmodel = load('./data/pnlssout_try0.mat', 'model');

Xpowers = pnlssmodel.model.xpowers(:,1:2);

[linmodel.Phi, linmodel.D] = eig(linmodel.model.A);
linmodel.E = zeros(2, size(Xpowers,1));
linmodel.E(2, find(Xpowers(:,1)==3 & Xpowers(:,2)==0)) = 1; % linmodel.model.nlcof.coef;
linmodel.Tmat = NLCOEF_TFMMATS(linmodel.Phi, Xpowers);

linmodel.phimodel = struct('A', linmodel.D, 'B', inv(linmodel.Phi)*linmodel.model.B, ...
    'C', linmodel.model.C*linmodel.Phi, 'D', linmodel.model.D, ...
    'E', inv(linmodel.Phi)*linmodel.E*linmodel.Tmat);


[pnlssmodel.Phi, pnlssmodel.D] = eig(logm(pnlssmodel.model.A)*fs);
pnlssmodel.Tmat = NLCOEF_TFMMATS(pnlssmodel.Phi, Xpowers);
Btmp = inv(pnlssmodel.model.A-eye(2))*logm(pnlssmodel.model.A)*fs*pnlssmodel.model.B;

pnlssmodel.phimodel = struct('A', pnlssmodel.D, ...
    'B', inv(pnlssmodel.Phi)*Btmp, ...
    'C', pnlssmodel.model.C*pnlssmodel.Phi, 'D', pnlssmodel.model.D, ...
    'E', inv(pnlssmodel.Phi)*pnlssmodel.model.E*pnlssmodel.Tmat);

El = linmodel.phimodel.E;
Ep = pnlssmodel.phimodel.E;