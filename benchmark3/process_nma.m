function nma=load_nma(nma, Ntd)

% dimensionless time - exclude last point, ie. tau = [0, 2pi[
tau = linspace(0,2*pi,Ntd+1); tau = tau(1:end-1);

PHI_L2 = nma.PHI_L2;
n = nma.n;
Psi = nma.X(1:end-3,:);
Om = nma.X(end-2,:);
del_HB = nma.X(end-1,:);
log10a = nma.X(end,:);
a = 10.^log10a;
Q_HB = Psi.*repmat(a,size(Psi,1),1);

w_L2 = zeros(Ntd,length(Om)); 
a_q = zeros(n,length(Om));
arms = 0;
% loop over all modes
for k = 1:n
    % get complex coefficients
    Qc = [Q_HB(k,:);Q_HB(n+k:2*n:end,:)-1i*Q_HB(2*n+k:2*n:end,:)];
    w_mode = PHI_L2(k)*real(exp(1i*tau(:)*(0:nma.H))*Qc);
    w_L2= w_L2 + w_mode;
    
    q = real(exp(1i*tau(:)*(0:nma.H))*Qc);
    a_q(k,:) = (max((q))-min((q)))/2;
    
    Q_rms = sqrt(sum(Q_HB(k:n:end,:).^2))/sqrt(2);
    arms = arms + PHI_L2(k)*Q_rms;
end

nma.arms = arms;
nma.apeak = (max((w_L2))-min((w_L2)))/2;
nma.tau = tau;
nma.Ntd = Ntd;
nma.Om = Om;
nma.Qsc = Q_HB;
nma.Qc = Qc;
nma.a_q = a_q;
end