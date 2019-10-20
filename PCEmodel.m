vmax = [30 27.5];
vc = 25;
L = [5 30]; %L1 was 6, L2 was 18
kjam = 1/L(1);
kc = kjam/6;
w = (vc*pc)/(pjam-pc);

% fundamental relation
vl = @(xT,n) vmax(n) - ((vmax(n)-vc)/kc).*xT; %xT:Total density
vr = @(xT) w.*(kjam./xT-1);
dvl = @(xT,n,eta) ((vmax(n)-vc)/kc).*xT.*eta;
dvr = @(xT,eta) w.*(-kjam./xT.^2).*eta;
ql = @(x,xT,n) x.*v(xT,n);
qr = @(x,xT,n) x.*v(xT,n);

% Initial density
ux0 = zeros(length(x),2);
% pce
pce = [1 2];
% Calculate density
[U,U1,U2,t] = NLLF3(x,T,ux0,vl,dvl,q,pce);
