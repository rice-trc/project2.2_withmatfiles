% simulate the system using ode45.
%
% This requires that the odesys is written in such a way that it can handle
% adaptive time steps. This is not the case if the multisine is calculated
% a priori, fx. using PNLSS

close all
clearvars

%% Define system
% Fundamental parameters
Dmod = [.38 .12 .09 .08 .08]*.01;
Nmod = 5;
setup = 'New_Design_Steel';
thickness = .001;
[L,rho,E,om,PHI,~,gam] = beams_for_everyone(setup,Nmod,thickness);
PHI_L2 = PHI(L/2);

% Properties of the underlying linear system
M = eye(Nmod);
D = diag(2*Dmod(:).*om(:));
K = diag(om.^2);

% load nonlinear coefficients (can be found e.g. analytically)
fname = ['beam_New_Design_Steel_analytical_5t_' ...
    num2str(thickness*1000) 'mm.mat'];
[p, E] = nlcoeff(fname, Nmod);

% Fundamental harmonic of external forcing
Fex1 = gam;

nz = size(p,1);
n = Nmod;

%% multisine, using time domain formulation

R = 4;           % Realizations. (one for validation and one for testing)
P = 8;           % Periods, we need to ensure steady state
f1 = 5;          % low freq
f2 = 400;        % high freq
fs = 1500;       % 5*f2. Must be fs>2*f2. Nyquist freq, you know:)
f0 = 1;          % freq resolution. 
N = f2/f0;       % freq points
har = ones(N,1); % full multisine, eg. excite all lines.
A = 15;           % amplitude
lines = 2:N;     % excited lines

upsamp = 1;         % upsampling factor to ensure integration accuracy.
fsint = fs*upsamp;  % integration sampling frequency.
Pfilter = 1;        % extra period to avoid edge effects during low-pass filtering
P = P + Pfilter;

q0 = zeros(n,1);
u0 = zeros(n,1);
t1 = 0;
t2 = P/f0;
Nt = fsint/f0;             % time points per period
t = t1:1/fsint:t2-1/fsint; % time vector. ode45 interpolates output
freq = (0:Nt-1)/Nt*fsint;  % frequency content
nt = length(t);

u = zeros(Nt,P,R);
y = zeros(Nt,P,R,n);
ydot = zeros(Nt,P,R,n);
tic
for r=1:R
    % predictable pseudo-random numbers
    rng(r);
    % multisine in time domain (sum of sines)
    phase = 2*pi*rand(N,1);
    fex = @(t) har'*A*cos(2*pi*(1:N)'*f0*t + phase) / sqrt(sum(har));

    par = struct('M',M,'C',D,'K',K,'p',p,'E',E,'fex',fex, 'amp', Fex1);
    [tout,Y] = ode45(@(t,y) sys(t,y, par), t,[q0;u0]);
 
    u(:,:,r) = reshape(fex(tout'), [Nt,P]);
    y(:,:,r,:) = reshape(Y(:,1:n), [Nt,P,n]);
    ydot(:,:,r,:) = reshape(Y(:,n+1:end), [Nt,P,n]);
end
disp(['ode45 with multisine in time domain required ' num2str(toc) ' s.']);

save('ode45_multisine.mat','u','y','ydot','f1','f2','fs','f0','freq',...
    't','A','PHI_L2','lines')

%% show time series
r = 1;
Y = PHI_L2*reshape(y(:,:,r,:),[],n)';

figure;
plot(tout, Y,'k-')
% indicate periods
h1 = vline(t((1:r*P)*Nt),'--g');
% indicate realizations
h2 = vline(t((1:r)*Nt*P),'--k');set([h1 h2],'LineWidth',0.5)
xlabel('time (s)')
ylabel('magnitude')
title(['Multisine: ' num2str(r) ' realizations of ' num2str(P) ' periods of ' num2str(Nt) ' samples per period'])
% export_fig('fig/multisine_sim_time.pdf')

% show periodicity
Y1 = reshape(Y,[Nt,P,r]);
figure;
per = (Y1(:,1:end-1,1)-Y1(:,end,1)) / rms(Y1(:,1,1));
plot(t(1:Nt*(P-1)),db(per(:)),'k-')
% indicate periods
h1 = vline(t((1:r*(P-1))*Nt),'--g');
% indicate realizations
h2 = vline(t((1:r)*Nt*(P-1)),'--k');set([h1 h2],'LineWidth',0.5)
xlabel('time (s)')
ylabel('Relative error to last period (dB)')
title([num2str(Nt) ' samples per period'])
% export_fig('fig/multisine_periodicity_lowfreq.pdf')

% force signal for one period
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
% export_fig('fig/multisine_freq.pdf')


%% Low-pass filtering and downsampling.
if upsamp > 1
    drate = factor(upsamp);        % prime factor decomposition.
    for k=1:length(drate),
        Y = decimate(Y,drate(k),'fir');
    end %k
    u = downsample(u,upsamp);
    
    % Removal of the last simulated period to eliminate the edge effects
    % due to the low-pass filter.
    Y = Y(1:(P-1)*Nt);
    u = u(1:(P-1)*Nt);
    P = P-1;
end

