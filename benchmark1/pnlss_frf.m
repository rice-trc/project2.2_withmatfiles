clc
clear all
addpath('../src/mhbm_res/')

load('./data/ode45_multisine.mat', 'fs', 't', 'u')
load('./data/pnlssout_try0.mat','model');

%%
Nh = 3;
Ntd = 2^10;
Ur = zeros(2*Nh+1, 1);
Ur(2) = 1;

Uc = Ur([1 2:2:end]) - 1j*[0; Ur(3:2:end)];

%% Gradients Check
rng(1)
X = rand((2*Nh+1)*size(model.A,1),1);
Om = 10;
rng(toc)
hv = zeros(size([X;Om]));  hm = 1e-6; hi = randi(size([X;Om],1))
hv(hi) = 1;
[R, J] = mhbm_aft_residual_pnlss_discrete([X; Om], model.A, model.B, model.E, model.xpowers, 1/fs, Uc, Nh, Ntd);
Rp = mhbm_aft_residual_pnlss_discrete([X; Om]+hv*hm, model.A, model.B, model.E, model.xpowers, 1/fs, Uc, Nh, Ntd);
Rm = mhbm_aft_residual_pnlss_discrete([X; Om]-hv*hm, model.A, model.B, model.E, model.xpowers, 1/fs, Uc, Nh, Ntd);
[(Rp-Rm)/(2*hm) J(:,hi)]

%% Run with NLvib
Ws = 100*2*pi;
We = 600*2*pi;
ds = 10*2*pi;
Sopt = struct('jac', 'full', 'stepmax', 1e4);

X0 = rand((2*Nh+1)*size(model.A,1),1);

Xr = solve_and_continue(X0, @(X) mhbm_aft_residual_pnlss_discrete(X, model.A, model.B, model.E, model.xpowers, 1/fs, Uc, Nh, Ntd),...
    Ws, We, ds, Sopt);
Xc = Xr([1 2:2:end],:) - 1j*[Xr(1,:)*0; Xr(3:2:end, :)];
Sol = [];
model.F = model.E(1,:)*0;
Ws = Xr(end,:);
As = zeros(size(Xr(end,:)));
for i=1:size(Xr,2)
    Sol = mhbm_post_amplitude_pnlss(Sol, Xc(:,i), 1, Uc, model.C, model.D, model.F, model.xpowers, Nh, Ntd);
    As(i) = Sol.A;
end

Yr = kron(eye(2*Nh+1), model.C)*Xr(1:end-1,:);
YAmp = sqrt([1 0.5*ones(1,2*Nh)]*Yr.^2)

plot(Ws/2/pi, As, '.-')