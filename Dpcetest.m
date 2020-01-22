%Dynamic PCE test 2
vmax = [30 27.5];
vc = 25;
L = [6 18]; %L1 was 6, L2 was 18
kjam = 1/L(1);
kc = kjam/6;
w = (vc*kc)/(kjam-kc);
TH = [1 1.5];

% fundamental relation
% uv and udv are universial variables
v = @(xT,n) uv(xT,n,vmax,vc,L(1));
dxT = @(x1,x2,xT,n) udxT(x1,x2,xT,n,L,TH,vmax,vc,kjam,kc);
dv = @(x1,x2,xT,m,n) udv(xT,m,vmax,vc,L(1)).*dxT(x1,x2,xT,n);
q = @(x,xT,n) x.*v(xT,n);
xT = @(x1,x2) uxT(x1,x2,L,TH,vmax,vc,kjam,kc);
pce = @(xT,n) (L(n)+TH(n).*v(xT,n))./(L(1)+TH(1).*v(xT,1));

% Discretization
x0 = 0;
xend = 5;
dx = 0.1;
x = x0:dx:xend;
T = 288;%2880; %1200

% Initial density
pt = [0.1 0.5*kc;-0.1 kjam;4 kjam;4 0];
dis = [0.8 0.2];
ux0 = Up(x,pt)'*dis./[pce(Up(x,pt),1);pce(Up(x,pt),2)]';

%% Calculate density
%[U,U1,U2,t] = CU3(x,T,ux0,v,dv,q,xT);
[UCUT,U1CUT,U2CUT,tCUT] = CU4(x,T,ux0,v,dv,q,xT);
[ULFT,U1LFT,U2LFT,tLFT] = NLLF4(x,T,ux0,v,dv,q,xT);
% graph
figure()
imagesc(tLFT,x,ULFT)
colorbar()
set(gca, 'XLim', tLFT([1 end]), 'YLim', x([1 end]), 'YDir', 'normal')
xticks([600 800 1200])

testpt = floor(length(x)/2);
figure()
av = (q(U1(testpt,:),U(testpt,:),1)+q(U2(testpt,:),U(testpt,:),2))./U(testpt,:);
plot(U(testpt,:),av,'.')
xlabel('Total Density')
ylabel('Average speed')

figure()
qT = (q(U1(testpt,:),U(testpt,:),1)+q(U2(testpt,:),U(testpt,:),2));
plot(U(testpt,:),qT,'.')
xlabel('Density')
ylabel('Flow')

figure()
for i = 1:7
    hold on
    plot(x,U(:,i*2e3));
    legend();
end
