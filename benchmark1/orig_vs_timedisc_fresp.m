clc
clear all
addpath('../src/mhbm_res/')
addpath('../src/matlab/')
addpath('../src/nlvib/SRC/')
addpath('../src/nlvib/SRC/MechanicalSystems/')

%% Define system

% Fundamental parameters
Dmod = [.38 .12 .09 .08 .08]*.01;
Nmod = 1;
setup = 'New_Design_Steel';
thickness = .001;
[L,rho,E,om,PHI,~,gam] = beams_for_everyone(setup,Nmod,thickness);
PHI_L2 = PHI(L/2);

% load nonlinear coefficients (can be found e.g. analytically)
fname = ['beam_New_Design_Steel_analytical_5t_' ...
    num2str(thickness*1000) 'mm.mat'];
[p, E] = nlcoeff(fname, Nmod);
E = 0;

% Properties of the underlying linear system
M = eye(Nmod);
D = diag(2*Dmod(1:Nmod).*om(1:Nmod));
K = diag(om.^2);

% Fundamental harmonic of external forcing
Fex1 = gam;

% Define oscillator as system with polynomial stiffness nonlinearities
oscillator = System_with_PolynomialStiffnessNonlinearity(M,D,K,p,E,Fex1);

% Number of degrees of freedom
n = oscillator.n;

%% Frequency response of original system using HBM
analysis = 'FRF';
H = 7;              % harmonic order
N=2*3*H+1;
Ntd = 1e3;

% Analysis parameters
Om_s = 200*2*pi;      % start frequency
Om_e = 700*2*pi;     % end frequency

% Excitation levels
exc_lev = [10 40 60 80 100];
X = cell(size(exc_lev));
Sols = cell(size(exc_lev));
for iex=1:length(exc_lev)
    % Set excitation level
    oscillator.Fex1 = Fex1*exc_lev(iex);
    
    % Initial guess (solution of underlying linear system)
    Q1 = (-Om_s^2*M + 1i*Om_s*D + K)\oscillator.Fex1;
    y0 = zeros((2*H+1)*length(Q1),1);
    y0(length(Q1)+(1:2*length(Q1))) = [real(Q1);-imag(Q1)];
    qscl = max(abs((-om(1)^2*M + 1i*om(1)*D + K)\oscillator.Fex1));
    
    % Solve and continue w.r.t. Om
    ds = 50; % -> better for exc_lev = 50
        
    Dscale = [1e-6*ones(length(y0),1);(Om_s+Om_e)/2];
    Sopt = struct('Dscale',Dscale,'dynamicDscale',1,'jac','full','stepmax',1e4);
    X{iex} = solve_and_continue(y0,...
        @(X) HB_residual(X,oscillator,H,N,analysis),...
        Om_s,Om_e,ds,Sopt);
	Sols{iex} = [X{iex}(end,:)' PHI_L2*sqrt([1 0.5*ones(1,2*H)]*X{iex}(1:end-1,:).^2)'...
        atan2d(-X{iex}(3,:),X{iex}(2,:))'];
end

%% Compute frequency response of PNLSS identified model
load('./data/ode45_multisine_A25.mat', 'fs');
load('./data/pnlssout_A25.mat','model');
fs = 2^18;

% Continuous time model
ctmodel.A = [0 1;
            -M\K -M\D];
ctmodel.B = [0; inv(M)*Fex1];
ctmodel.C = [PHI_L2 0];
ctmodel.D = [0];
ctmodel.xpowers = [3 0 0];
ctmodel.E = [0; -M\E];

ctmodel.xpowers = model.xpowers;
ctmodel.E = [[0; -M\E] zeros(2, 9)];

% Discrete time model
dtmodel = ctmodel;
dtmodel.A = ctmodel.A/fs+eye(2)
dtmodel.B = ctmodel.B/fs;
dtmodel.E = ctmodel.E/fs;

dtmodel.A = expm(ctmodel.A/fs);
dtmodel.B = inv(ctmodel.A)*(dtmodel.A-eye(2))*ctmodel.B;
dtmodel.E = ctmodel.E/fs;

Ndpnlss = size(dtmodel.A,1);

% Forcing vector
Uc = zeros(H+1,1);
Uc(2) = 1;

ds = 1*2*pi;
dsmin = 0.001*2*pi;
dsmax = 50*2*pi;

Xpnlss = cell(length(exc_lev),1);
Solspnlss = cell(length(exc_lev),1);
for iex=1:length(exc_lev)
    Ff = exc_lev(iex);

    Xc = (exp(1i*Om_s/fs)*eye(size(dtmodel.A))-dtmodel.A)\(dtmodel.B*Ff);             % linear solution
    % Xc = ((1i*Om_s/fs)*eye(size(dtmodel.A))-dtmodel.A)\(dtmodel.B*Ff);             % linear solution
    X0 = [zeros(length(dtmodel.A),1);real(Xc);-imag(Xc);....
            zeros(2*(H-1)*length(dtmodel.A),1)];                  % initial guess
    
    Dscale = [mean(abs(Xc))*ones(length(X0),1);Om_s];
    Sopt = struct('ds',ds,'dsmin',dsmin,'dsmax',dsmax,'flag',1,'stepadapt',1, ...
            'predictor','tangent','parametrization','arc_length', ...
            'Dscale',Dscale,'jac','full', 'dynamicDscale', 1);

    fun_residual = ...
            @(XX) mhbm_aft_residual_pnlss_discrete(XX, dtmodel.A, dtmodel.B, dtmodel.E, dtmodel.xpowers, 1/fs, Uc*Ff, H, Ntd);
    Cfun_postprocess = {@(varargin) ...
            mhbm_post_amplitude_pnlss(varargin{:},Uc*Ff,dtmodel.C,dtmodel.D,zeros(1,length(dtmodel.E)),dtmodel.xpowers,H,Ntd)};
    fun_postprocess = @(Y) mhbm_postprocess(Y,fun_residual,Cfun_postprocess);

    [Xpnlss{iex},~,Sol] = solve_and_continue(X0, fun_residual,...
        Om_s, Om_e, ds, Sopt, fun_postprocess);
    Solspnlss{iex} = [Xpnlss{iex}(end,:)' [Sol.Apv]' [Sol.Aph1]'];
end


%% Plot
figure(1)
clf()

figure(2)
clf()
colos = distinguishable_colors(length(exc_lev));
aa = gobjects(size(exc_lev));
for iex=1:length(exc_lev)
    figure(1)
    plot(Sols{iex}(:,1)/2/pi, Sols{iex}(:,2), '-', 'Color', colos(iex,:)); hold on
    plot(Solspnlss{iex}(:,1)/2/pi, Solspnlss{iex}(:,2), '.--', 'Color', colos(iex,:))
    
    figure(2)
    aa(iex) = plot(Sols{iex}(:,1)/2/pi, Sols{iex}(:,3), '-', 'Color', colos(iex,:)); hold on
    plot(Solspnlss{iex}(:,1)/2/pi, Solspnlss{iex}(:,3), '.--', 'Color', colos(iex,:))
    legend(aa(iex), sprintf('F = %.2f', exc_lev(iex)));
end

figure(1)
set(gca,'Yscale', 'log')
xlim(sort([Om_s Om_e])/2/pi)
xlabel('Forcing frequency \omega (Hz)')
ylabel('RMS response amplitude (m)')
% savefig(sprintf('./fig/pnlssfrf_A%d_Amp.fig',Alevel))
% print(sprintf('./fig/pnlssfrf_A%d_Amp.eps',Alevel), '-depsc')

figure(2)
xlim(sort([Om_s Om_e])/2/pi)
xlabel('Forcing frequency \omega (Hz)')
ylabel('Response phase (degs)')
legend(aa(1:end), 'Location', 'northeast')
% savefig(sprintf('./fig/pnlssfrf_A%d_Phase.fig',Alevel))
% print(sprintf('./fig/pnlssfrf_A%d_Phase.eps',Alevel), '-depsc')