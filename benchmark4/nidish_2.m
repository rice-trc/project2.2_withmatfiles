clc
clear all

addpath('./NLvib_v1.1/SRC/')
addpath('./NLvib_v1.1/SRC/MechanicalSystems/')

%% System definition
len = 0.70;
hgt = 0.03;
thk = hgt;
E   = 185e9;
rho = 7830.0;
BCs = 'clamped-free';

Nn = 8;
beam = FE_EulerBernoulliBeam(len, hgt, thk, E, rho, BCs, Nn);
Fex1 = zeros(size(beam.M,1),1);  Fex1(end-1) = 1;
% Nonlinearity
Nnl = 4;
dir = 'trans';
kt  = 1.3e6;
muN = 1;
add_nonlinear_attachment(beam, Nnl, dir, 'elasticdryfriction', ...
    'stiffness', kt, 'friction_limit_force', muN, ...
    'ishysteretic', true);

%% Linearized limit cases
% Slipped
[Vsl, Dsl] = eig(beam.K, beam.M);
[Dsl, si] = sort(sqrt(diag(Dsl)));
Vsl = Vsl(:, si);  Vsl = Vsl./sqrt(diag(Vsl'*beam.M*Vsl)');

% Stuck
Knl = zeros(size(beam.M));
nli = find(beam.nonlinear_elements{1}.force_direction);
Knl(nli, nli) = kt;
Kst = beam.K + Knl;
[Vst, Dst] = eig(Kst, beam.M);
[Dst, si] = sort(sqrt(diag(Dst)));
Vst = Vst(:, si); Vst = Vst./sqrt(diag(Vst'*beam.M*Vst));

%% Rayleigh damping
% Desired
zs = [8e-3; 2e-3];
ab = [ones(length(zs),1) Dst(1:length(zs)).^2]\(2*zs.*Dst(1:length(zs)));

beam.D = ab(1)*beam.M + ab(2)*Kst;

Zetas = diag(Vst'*beam.D*Vst)./(2*Dst);

%% Frequency Response
analysis = 'FRF';
Nh = 9;
Nhc = 2*Nh+1;
Nd = size(beam.M, 1);
Nt = 2^10;

Ws = 200;
We = 500;
ds = (We-Ws)/100;
dsmax = (We-Ws)/50;

Dscale = [1e-4*ones(Nd*Nhc,1); (Dst(1)+Dsl(1))/2];

Fas = [0.01 0.05 0.1 0.2 0.3 0.6];
Xcont = cell(size(Fas));
for k=1:length(Fas)
    beam.Fex1(end-1) = Fas(k);  % Forcing from last node
    H1tmp = ((Kst-Ws^2*beam.M)+1j*(Ws*beam.D))\beam.Fex1;
    X0 = zeros(Nd*Nhc, 1); X0(Nd+(1:2*Nd)) = [real(H1tmp); -imag(H1tmp)];
    
    Sopt = struct('jac','full','stepmax',1000,'MaxFfunEvals',500, ...
        'dsmax', dsmax, 'Dscale', Dscale, 'dynamicDscale', 1); %, ...
%         'parametrization', 'normal');
    Xcont{k} = solve_and_continue(X0, ...
        @(X) HB_residual(X, beam, Nh, Nt, analysis), Ws, We, ds, Sopt);
end

%% EPMC NMA
analysis = 'nma';
imod = 1;
log10a_s = -7;
log10a_e = -3.8;
inorm = Nd-1;

X0 = zeros(Nd*Nhc+2, 1);
X0(Nd+(1:Nd)) = Vst(:,imod);
X0(end-1) = Dst(imod);
X0(end) = 2*Zetas(imod)*Dst(imod);

Dscale = [1e-1*ones(Nd*Nhc,1); Dst(1,1); 1.0; 1.0];

beam.Fex1 = beam.Fex1*0;

ds = 0.01;
dsmax = 0.05;
Sopt = struct('jac', 'full', 'Dscale', Dscale, 'dynamicDscale', 1);
% Sopt = struct('jac','full','stepmax',1000,'MaxFfunEvals',500, ...
% 	'dsmax', dsmax); %, ...
%         'parametrization', 'normal');
fscl = mean(abs(beam.K*Vst(:,imod)));
Xbb = solve_and_continue(X0, ...
    @(X) HB_residual(X, beam, Nh, Nt, analysis, inorm, fscl), ...
    log10a_s, log10a_e, ds, Sopt);

%% Plotting with NM-ROM
figure(1);
clf();
colos = distinguishable_colors(length(Fas), 'k');

aa = gobjects(length(Fas)+1,1);
aa(1) = semilogy(Xbb(end-2,:)/(2*pi), 10.^Xbb(end,:).*sqrt([1 0.5*ones(1,Nhc-1)]*Xbb((Nd-1):Nd:end-2,:).^2), ...
        'k+-', 'LineWidth', 2); hold on
legend(aa(1), 'EPMC')

a_NMA = 10.^(Xbb(end,:)');
p2 = (Xbb(end-2,:).^2-2*(Xbb(end-1,:).*Xbb(end-2,:)).^2)';
om4 = (Xbb(end-2,:).^4)';
Phi_HB = Xbb(Nd+(1:Nd),:)-1j*Xbb(2*Nd+(1:Nd),:);
Fsc = (abs(Phi_HB'*Fex1)./a_NMA).^2;
mAmps = sqrt(0.5)*a_NMA.*abs(Phi_HB(Nd-1,:))';

livs = 10;
for k=1:length(Fas)
    om1 = sqrt(p2 + sqrt(p2.^2-om4+Fsc*Fas(k)^2));
    om2 = sqrt(p2 - sqrt(p2.^2-om4+Fsc*Fas(k)^2));

    ris1 = find(imag(om1)==0);  % real indices
    ris2 = find(imag(om2)==0);  % real indices
    semilogy(om1(ris1)/(2*pi), mAmps(ris1), '--', 'Color', colos(k,:)); hold on
    semilogy(om2(ris2)/(2*pi), mAmps(ris2), '--', 'Color', colos(k,:)); hold on
    
    ivs = fix(linspace(1, length(ris1), livs));
	semilogy(om1(ris1(ivs))/(2*pi), mAmps(ris1(ivs)), ...
        'o', 'Color', colos(k,:), 'MarkerFaceColor', colos(k,:)); hold on
    ivs = fix(linspace(1, length(ris2), livs));
    semilogy(om2(ris2(ivs))/(2*pi), mAmps(ris2(ivs)), 'o', ...
        'Color', colos(k,:), 'MarkerFaceColor', colos(k,:)); hold on

    aa(k+1) = semilogy(Xcont{k}(end,:)/(2*pi), sqrt([1 0.5*ones(1,Nhc-1)]*Xcont{k}((Nd-1):Nd:end-1,:).^2), ...
        '-', 'LineWidth', 2, 'Color', colos(k,:)); hold on
    legend(aa(k+1), sprintf('F = %.2f N', Fas(k)));
end
legend(aa(1:end), 'Location', 'northeast')
xlim([30 85])
xlabel('Frequency (Hz)')
ylabel('Response Amplitude (m)')
% print('./FIGS/NMROMRES.eps', '-depsc')

xlim([58 68])
ylim([3e-7 1e-5])
legend('hide')
% print('./FIGS/NMROMRES_ZOOM1.eps', '-depsc')

xlim([41 59])
ylim([5e-7 1e-4])
% print('./FIGS/NMROMRES_ZOOM2.eps', '-depsc')