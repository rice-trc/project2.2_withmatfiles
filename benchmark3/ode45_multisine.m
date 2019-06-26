% simulate the system using ode45.
%
% This requires that the odesys is written in such a way that it can handle
% adaptive time steps. This is not the case if the multisine is calculated
% a priori, fx. using PNLSS

close all
clear all

srcdir = '../src/matlab';
addpath(genpath(srcdir));

dataname = 'ms_odd_random';
savedata = true;
savefig = true;

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
modelname  = ['data/','beam_msh_80_4_1_3t_steel_',...
    num2str(thickness*1000),'mm_R',num2str(R),'m.mat'];
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

%% multisine, using time domain formulation

R = 1;           % Realizations. (one for validation and one for testing)
P = 5;           % Periods, we need to ensure steady state
f1 = 5;          % low freq
f2 = 450;        % high freq
fs = 1500;       % 5*f2. Must be fs>2*f2. Nyquist freq, you know:)
f0 = 1;          % freq resolution. 
N = f2/f0;       % freq points
A = 15;          % amplitude
Nt = fs/f0;      % time points per period
% set the type of multisine
ms_type = 'random_odd';  % can be 'full', 'odd' or 'random_odd'
ms = multisine_lines(N,Nt, ms_type);

q0 = zeros(n,1);
u0 = zeros(n,1);
t1 = 0;
t2 = P/f0;

t = t1:1/fs:t2-1/fs;   % time vector. ode45 interpolates output
freq = (0:Nt-1)/Nt*fs; % frequency content
nt = length(t);        % total number of points

u = zeros(Nt,P,R);
y = zeros(Nt,P,R,n);
ydot = zeros(Nt,P,R,n);
phases = zeros(N,R);

tic
for r=1:R
    % predictable pseudo-random numbers
    rng(r);
    % multisine in time domain (sum of sines)
    phase = 2*pi*rand(N,1);
    fex = @(t) ms.har'*A*cos(2*pi*(1:N)'*f0*t + phase) / sqrt(sum(ms.har));

    par = struct('M',M,'C',D,'K',K,'p',p,'E',E,'fex',fex, 'amp', Fex1);
    [tout,Y] = ode45(@(t,y) sys(t,y, par), t,[q0;u0]);
 
    u(:,:,r) = reshape(fex(tout'), [Nt,P]);
    y(:,:,r,:) = reshape(Y(:,1:n), [Nt,P,n]);
    ydot(:,:,r,:) = reshape(Y(:,n+1:end), [Nt,P,n]);
    phases(:,r) = phase;
end
disp(['ode45 with multisine in time domain required ' num2str(toc) ' s.']);

if savedata
save(['data/',dataname],'u','y','ydot','f1','f2','fs','f0','freq',...
    't','A','PHI_L2','ms','phases')
end

%% init plotting
hfig = {};
% select realization to do the plotting for
r = 1;

%% show time series
yres = PHI_L2*reshape(y(:,:,r,:),[],n)';

figure;
plot(tout, yres,'k-')
% indicate periods
h1 = vline(t((1:r*P)*Nt),'--g');
% indicate realizations
h2 = vline(t((1:r)*Nt*P),'--k');set([h1 h2],'LineWidth',0.5)
xlabel('time (s)')
ylabel('magnitude')
title(['Multisine: ' num2str(r) ' realizations of ' num2str(P),...
    ' periods of ' num2str(Nt) ' samples per period'])

hfig(end+1) = {{gcf, "ms_time"}};


%% show periodicity

figure;
per = (y(:,1:end-1,r)-y(:,end,r)) / rms(y(:,1,r));
plot(t(1:Nt*(P-1)),db(per(:)),'k-')
% indicate periods
h1 = vline(t((1:r*(P-1))*Nt),'--g');
% indicate realizations
h2 = vline(t((1:r)*Nt*(P-1)),'--k');set([h1 h2],'LineWidth',0.5)
xlabel('time (s)')
ylabel('Relative error to last period (dB)')
title([num2str(Nt) ' samples per period'])

hfig(end+1) = {{gcf, "ms_periodicity"}};


%% show force signal for one period
% plot only half spectrum(single sided spectrum)
dft = fft(u(:,1,r));
nt = length(dft)/2+1;
xdft = dft(1:nt);

figure; subplot(2,1,1)
plot(freq(1:nt),db(xdft),'-*')
xlabel('frequency (Hz)')
ylabel('magnitude (dB)')
title('FFT of one period of the multisine realizations')

subplot(2,1,2)
plot(freq(1:nt),angle(xdft),'-*')
xlabel('frequency (Hz)')
ylabel('phase (rad)')
title('FFT of one period of the multisine realizations')

hfig(end+1) = {{gcf, "ms_exfreq"}};


%% show detect type of nonlinearity

% output signal for one period
% plot only half spectrum(single sided spectrum)
dft = fft(y(:,end,r));
nt = length(dft)/2+1;
xdft = dft(1:nt);
Ydot = fft(ydot(:,end,r));

figure;
hold on
plot(freq(ms.lines),db(xdft(ms.lines)),'.k')
scatter(freq(ms.non_odd),db(xdft(ms.non_odd)),[])
scatter(freq(ms.non_even),db(xdft(ms.non_even)),[])
hold off

xlabel('frequency (Hz)')
ylabel('magnitude (dB)')
title('FFT of one period of the output spectrum')
legend('excited','non\_odd','non\_even')

hfig(end+1) = {{gcf, "ms_nl_detection"}};

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