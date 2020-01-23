%Constant PCE
vmax = [30 27.5];
vc = 25;
L = [6 18]; %L1 was 6, L2 was 18
kjam = 1/L(1);
kc = kjam/6;
w = (vc*kc)/(kjam-kc);
pce = [1 L(2)/(L(1))];

% fundamental relation
% uv and udv are universial variables
v = @(xT,n) uv(xT,n,vmax,vc,L(1));
dv = @(xT,m,n) udv(xT,m,vmax,vc,L(1))*pce(n);
q = @(x,xT,n) x.*v(xT,n);

% Discretization
x0 = -5500;
xend = 10000;%did 1000
dx = 5;
x = x0:dx:xend;
T = 1200; %1200

% Initial density
dis = [0.8 0.2]./pce;
pt = [-2000 0.5*kc;-2000 kjam;0 kjam;0 0];
ux0 = Up(x,pt)'*dis;

%% Calculate density
[UCUP,U1CUP,U2CU,tCUP] = CU3(x,T,ux0,v,dv,q,pce);
[ULFP,U1LF,U2LF,tLFP] = NLLF3(x,T,ux0,v,dv,q,pce);

% graph
figure()
imagesc(tLFP,x,ULFP)
colorbar()
set(gca, 'XLim', tLFP([1 end]), 'YLim', x([1 end]), 'YDir', 'normal')
xticks([600 800])

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
    plot(x,U(:,i*1e3));
    legend();
end
dxCU = 25;
xCU = x0:dxCU:xend;

dxLF = 5;
xLF = x0:dxLF:xend;

plot(xCU,UCUP(:,end),'b-')
hold on
plot(xLF,ULFP(:,end),'r-')

plot(xCU,UCU(:,end),'k*')
hold on
plot(xLF,ULF(:,end),'y*')
legend("CUP","LFP","CU","LF")