clear; clc;
% close all;
% srcpath = '../src/nlvib';
% rmpath(genpath(srcpath));
% addpath(genpath(srcpath));

addpath('00_SRC');
addpath('00_SRC/MechanicalSystems');

srcpath = '../src/matlab';
addpath(genpath(srcpath));

savedata = true;
savename = '';

%% Define system

% Fundamental parameters
Dmod = [.38 .09 .08]*.01; % first, third, fifth bending modes of flat beam
Nmod = 3;
setup = 'New_Design_Steel';

thickness = .001;
R=3;

[L,rho,E,~,PHI,~,gam] = beams_for_everyone(setup,Nmod*2-1,thickness);
PHI_L2 = PHI(L/2);
PHI_L2 = PHI_L2(1:2:end);
gam = gam(1:2:end);

% load nonlinear coefficients (can be found e.g. via IC-method)
modelname  = ['beam_msh_80_4_1_3t_steel_' num2str(thickness*1000) 'mm_R' num2str(R) 'm.mat'];
[p, E] = nlcoeff(modelname, Nmod);

% also om is needed for linear model, so we load the model again
load(modelname);
om = model.omega;

% Properties of the underlying linear system
M = eye(Nmod);
D = diag(2*Dmod(:).*om(:));
K = diag(om.^2);


% Fundamental harmonic of external forcing
Fex1 = gam;

% Define oscillator as system with polynomial stiffness nonlinearities
oscillator = System_with_PolynomialStiffnessNonlinearity(M,D,K,p,E,Fex1);

% Number of degrees of freedom
n = oscillator.n;

%% NMA
analysis = 'NMA';
H=20;
N=2*3*H+1;

% Linear modal analysis
[PHI_lin,OM2] = eig(oscillator.K,oscillator.M);
[om_lin,ind] = sort(sqrt(diag(OM2)));
PHI_lin = PHI_lin(:,ind);

imod = 1;           % mode to be analyzed
log10a_s = -8;    % start vibration level (log10 of modal mass)
log10a_e = -3.2;       % end vibration level (log10 of modal mass)
inorm = 1;          % coordinate for phase normalization

% Initial guess vector y0 = [Psi;om;del], where del is the modal
% damping ratio, estimate from underlying linear system
om = om_lin(imod); phi = PHI_lin(:,imod);
Psi = zeros((2*H+1)*Nmod,1);
Psi(Nmod+(1:Nmod)) = phi;
x0 = [Psi;om;0];

ds      = .15;
Sopt    = struct('Dscale',[1e-6*ones(size(x0,1)-2,1);1;1e-1;1],...
    'dynamicDscale',1,'stepmax',5e4);
% Sopt    = struct('stepmax',5e4);
[X,Solinfo,Sol] = solve_and_continue(x0,...
    @(X) HB_residual(X,oscillator,H,N,analysis,inorm),...
    log10a_s,log10a_e,ds, Sopt);


nmadata = ws2struct('Solinfo','Sol','X','ds','H','n','imod','oscillator',...
    'Dmod','Nmod','om','gam','thickness','setup','modelname','R','PHI_L2');
if savedata
    save([savename,'nma.mat'],'-struct','nmadata')
end

%%
Xm1 = X(1:oscillator.n:end-3, :);
Xm1c = Xm1([1 2:2:end], :) - 1j*[Xm1(1,:)*0; Xm1(3:2:end,:)];
Xm1a = abs(Xm1c);

figure(5);
clf()
aa = gobjects(H, 1);
for nh=1:H
    aa(nh) = semilogy(X(end-2,:)/2/pi, Xm1a(1+nh,:), '-'); hold on
    legend(aa(nh), sprintf('%d', nh))
end
legend(aa(1:nh))
%%
% hb_plot