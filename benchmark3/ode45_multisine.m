% Simulate the system using ode5.
%
% To truly ensure periodicity(ie. steady state), it is important that you
% DO NOT use matlabs ode45 solver, as it is a variable step-size solver.
% A fixed step-size ode solver solver seeems to be working.


% close all
clear all

srcdir = '../src/matlab';
addpath(genpath(srcdir));

dataname = 'ms_full';
savedata = true;
savefig = true;

benchmark = 1;

%% Define system
switch benchmark
    case 1
        % Excitation levels. Note 50 is too high!
        exc_lev = [1,10,20,30];
        % low/high freq
        f1 = 50;
        f2 = 405;
    case 2
        exc_lev = [1,15,50];
        f1 = 50;
        f2 = 405;

    case 3
        exc_lev = [1,15,50];
        f1 = 50;
        f2 = 405;
end

sys = load_sys(benchmark);
n = sys.Nmod;

%% multisine, using time domain formulation

R  = 5;           % Realizations. (one for validation and one for testing)
P  = 10;           % Periods, we need to ensure steady state

% select Nt satisfying fs as 80% of Nyquist freq. Does not work yet.
% f0 = 0.1;
% N  = ceil((f2-f1)/f0);      % freq points
% Nt = ceil(1/0.8 * f2*2/f0); % Time points per cycle
% fs = Nt*f0;                 % Samping frequency - implicit given

N  = 2e3;         % freq points
f0 = (f2-f1)/N;
Nt = 2^14;      % Time per cycle
fs = Nt*f0;     % Samping frequency

% set the type of multisine
ms_type = 'full';  % can be 'full', 'odd' or 'random_odd'

if fs/2 <= f2
    error('Increase sampling rate! Current f0:%g, requires atleast Nt:%g',f0,f2*2/f0);
end

q0   = zeros(n,1);
u0   = zeros(n,1);
t1   = 0;
t2   = P/f0;
% time vector.
t    = linspace(t1, t2, Nt*P+1);  t(end) = [];
% t1 = t1:1/fs:t2-1/fs;
freq = (0:Nt-1)/Nt*fs;   % frequency content
nt = length(t);          % total number of points

fprintf('running ms benchmark:%d. R:%d, P:%d, Nt:%d, fs:%g, f0:%g\n',...
    benchmark,R,P,Nt,fs,f0);
for A = exc_lev
u = zeros(Nt,P,R);
y = zeros(Nt,P,R,n);
ydot = zeros(Nt,P,R,n);
MS = cell(R, 1);

tic
for r=1:R
    % multisine force signal
    [fex, MS{r}] = multisine(f1, f2, N, A, [], [], r);

    par = struct('M',sys.M,'C',sys.D,'K',sys.K,'p',sys.p,'E',sys.E,'fex',fex, 'amp', sys.Fex1);
%     [tout,Y] = ode45(@(t,y) odesys(t,y, par), t,[q0;u0]);
    Y = ode5(@(t,y) odesys(t,y, par), t,[q0;u0]);
 
    u(:,:,r) = reshape(fex(t'), [Nt,P]);
    y(:,:,r,:) = reshape(Y(:,1:n), [Nt,P,n]);
    ydot(:,:,r,:) = reshape(Y(:,n+1:end), [Nt,P,n]);
end
disp(['ode5 with multisine in time domain required ' num2str(toc) ' s.']);

if savedata
    save(sprintf('data/b%d_A%d_%s',benchmark,A,dataname),'u','y','ydot',...
        'f1','f2','fs','freq','t','A','sys','MS')
end

% only plot if it's supported (ie if we're not running it from cli)
% https://stackoverflow.com/a/30240946
if usejava('jvm') && feature('ShowFigureWindows')
    % phi: convert to physical coordinate
    phi = sys.PHI([sys.L/2]);
    ms_plot(t,y,u,freq,MS{r}.lines,phi,benchmark,A,dataname,savefig)
end

end
