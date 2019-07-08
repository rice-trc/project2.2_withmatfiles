clc
clear variables

addpath('../src/pnlss');

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

Om_s = 200*2*pi;      % start frequency
Om_e = 700*2*pi;     % end frequency


%% state space

fs = 2^16;

% Continuous time model
ctmodel.A = [0 1;
            -M\K -M\D];
ctmodel.B = [0; M\Fex1];
ctmodel.C = [PHI_L2 0];
ctmodel.D = [0];
% xpowers: each row gives the power of [x1, x2, u], which is then
% multiplied with corresponding coefficient in E.
% E ∈ [n, nx], xpowers ∈ [nx, n+1]
ctmodel.xpowers = [3 0 0];
ctmodel.E = [0; -M\E];
ctmodel.F = 0;
ctmodel.ypowers = [3 0 0];

% analytical discretization for A,B
dtmodel = ctmodel;
dtmodel.A = expm(ctmodel.A/fs);
dtmodel.B = ctmodel.A\(dtmodel.A-eye(2))*ctmodel.B;
dtmodel.E = ctmodel.E/fs;
sys = dtmodel;

% parameters for PNLSS. Initialize as linear. We change this later
% nx/ny is not needed for time simulation
nx = []; ny = []; T1 = 0; T2= 0;
model = fCreateNLSSmodel(sys.A,sys.B,sys.C,sys.D,nx,ny,T1,T2);
model.E = dtmodel.E;
model.F = dtmodel.F;
model.xpowers = dtmodel.xpowers;
model.ypowers = dtmodel.ypowers;
%% simualte model -- From here is new!
% direction to sweep
dir = 1; % or -1

% We see each increment in frequency as a new realization.
Om_vec = linspace(Om_s, Om_e, 50);  % Vector of excited frequencies
if dir == -1
    Om_vec = fliplr(Om_vec);
end
R = length(Om_vec);
P = 10;
Nt = 2^12;                          % points per period
t = zeros(length(Om_vec), Nt*P+1);  % end time depends on forcing freq
for i = 1:length(Om_vec)
    t(i,:) = linspace(0, P*2*pi/Om_vec(i), Nt*P+1);
end
t(:,end) = [];

% transient periods to simulate and remove to ensure steady state
NTrans = 10*Nt;  % Add 10 periods before the start
% Number of transient samples and starting indices of each realization
T1 = [NTrans 1+(0:Nt*P:(R-1)*Nt*P)]; 
model.T1 = T1;

% amplitudes to loop over
Avec = [10,40,80,100];
a_max = zeros(length(Om_vec),length(Avec));
j = 1;
for A = Avec
    u = A*sin(Om_vec(:).*t);
    u = reshape(u',[],1); % Nt*P*R x m
    
    y = fFilterNLSS(model,u); % Modeled output signal
    
    % calculate max amplitude
    yper = reshape(y, [Nt*P,R]);
    
    % take amplitude as max-min over last 5 periods
    nper = 5;
    for i = 1:length(Om_vec)
        ymax = max(yper(end-nper*Nt-1:end,i));
        ymin = min(yper(end-nper*Nt-1:end,i));
        a_max(i,j) = abs(0.5*(ymax-ymin));
    end
    j = j+1;
    fprintf('A: %g\n',A)
end

save(sprintf('data/stepped_sweep_dir%d',dir),'Om_vec','a_max')


%% plot max amplitudes -- poor mans frf
figure
hold on
for i = 1:length(Avec)
    plot(Om_vec/2/pi, a_max(:,i), 'x')
end
xlabel('Frequency (Hz)')
ylabel('Max amplitude (m)')

%% Plot time signal
y1 = reshape(y, [Nt*P,R]);

figure
hold on
for i = 1:length(Om_vec)
    plot(t(i,:), y1(:,i))
end
xlabel('Time (s)')
ylabel('Amplitude (m)')

%% plot freq content of last 5 periods

Y = fft(y1(end-Nt*5-1:end,:),[],1);  % along 1st dim
ndft = size(Y,1)/2;
figure
hold on
for i = 1:length(Om_vec)
plot(db(abs(Y(1:ndft,i))))
end

%% Check periodicity
yper = reshape(y, [Nt,P,R]);
r = 1;
figure;
per = (yper(:,1:end-1,r)-yper(:,end,r)) / rms(yper(:,end,r));
plot(t(r,1:Nt*(P-1)),db(per(:)),'k-')
% indicate periods
h1 = vline(t(r,(1:Nt:(P-1)*Nt)),'--g'); set([h1],'LineWidth',0.5)
xlabel('time (s)')
ylabel('Relative error to last period (dB)')
title([num2str(Nt) ' samples per period'])

fprintf('The mean error between last two periods: %g\n',mean(per(:,end)))

