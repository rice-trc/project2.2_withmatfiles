clc
clear all
addpath('../src/mhbm_res/')
addpath('../src/nlvib/SRC/')

load('./data/ode45_multisine.mat', 'fs', 't', 'u')
load('./data/pnlssout_try0.mat','model');

%%
Nd = size(model.A,1);
Nh = 3;
Ntd = 2^10;
Ur = zeros(2*Nh+1, 1);
Ur(2) = 1;

F = 5;

Uc = Ur([1 2:2:end]) - 1j*[0; Ur(3:2:end)];

Ur = Ur*F;
Uc = Uc*F;

% %% Gradients Check
% rng(1)
% X = rand((2*Nh+1)*size(model.A,1),1);
% Om = 10;
% rng(toc)
% hv = zeros(size([X;Om]));  hm = 1e-6; hi = randi(size([X;Om],1))
% hv(hi) = 1;
% [R, J] = mhbm_aft_residual_pnlss_discrete([X; Om], model.A, model.B, model.E, model.xpowers, 1/fs, Uc, Nh, Ntd);
% Rp = mhbm_aft_residual_pnlss_discrete([X; Om]+hv*hm, model.A, model.B, model.E, model.xpowers, 1/fs, Uc, Nh, Ntd);
% Rm = mhbm_aft_residual_pnlss_discrete([X; Om]-hv*hm, model.A, model.B, model.E, model.xpowers, 1/fs, Uc, Nh, Ntd);
% [(Rp-Rm)/(2*hm) J(:,hi)]

%% Run with NLvib
Ws = 320*2*pi;
We = 240*2*pi;
ds = 1;
% X0 = rand((2*Nh+1)*size(model.A,1),1);

Xc = (exp(1i*Ws/fs)*eye(size(model.A))-model.A)\(model.B*F);             % linear solution
X0 = [zeros(length(model.A),1);real(Xc);-imag(Xc);....
        zeros(2*(Nh-1)*length(model.A),1)];                  % initial guess
    
Sopt = struct('ds',ds,'dsmin',ds/5,'dsmax',5*ds,'flag',1,'stepadapt',1, ...
        'predictor','tangent','parametrization','arc_length', ...
        'Dscale',[10*ones(length(X0),1);Ws],'jac','full');

fun_residual = ...
        @(XX) mhbm_aft_residual_pnlss_discrete(XX, model.A, model.B, model.E, model.xpowers, 1/fs, Uc, Nh, Ntd);
Cfun_postprocess = {@(varargin) ...
        mhbm_post_amplitude_pnlss(varargin{:},Uc,model.C,model.D,zeros(1,length(model.E)),model.xpowers,Nh,Ntd)};
fun_postprocess = @(Y) mhbm_postprocess(Y,fun_residual,...
        Cfun_postprocess);

[Xr,~,Sol] = solve_and_continue(X0, fun_residual,...
    Ws, We, ds, Sopt, fun_postprocess);

%%
Sopt = struct('jac', 'full', 'stepmax', 1e4);
X0 = rand(Nd*(2*Nh+1),1);
Xr = solve_and_continue(X0, ...
    @(X) mhbm_aft_residual_pnlss_discrete(X, model.A, model.B, model.E, model.xpowers, 1/fs, Uc, Nh, Ntd),...
    Ws, We, ds, Sopt);
Xc = zeros(size(model.A,1)*(Nh+1), size(Xr,2));
Xc(1:Nd,:) = Xr(1:Nd,:);
for hi=1:Nh
    Xc(hi*Nd+(1:Nd),:) = Xr(Nd+(hi-1)*2*Nd+(1:Nd),:)-1j*Xr(2*Nd+(hi-1)*2*Nd+(1:Nd),:);
end
Sol = [];
Ws = Xr(end,:);
As = zeros(size(Xr(end,:)));
for i=1:size(Xr,2)
    Sol = mhbm_post_amplitude_pnlss(Sol, Xc(:,i), 1, Uc, model.C, model.D, model.F, model.ypowers, Nh, Ntd);
    As(i) = Sol.A;
end

As = [Sol.Amax];

Yr = kron(eye(2*Nh+1), model.C)*Xr(1:end-1,:);
YAmp = sqrt([1 0.5*ones(1,2*Nh)]*Yr.^2);

% testl = Xr(3,:) - 1j*Xr(4,:);
% plot(Xr(end,:)/2/pi, abs(testl), '.-')

hold on;
plot(Xr(end,:)/2/pi, As, '.-')