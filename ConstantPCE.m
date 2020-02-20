%Constant PCE
clear all
vmax = [30 27.5];%30 27.5
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
xT = @(x1,x2) x1*pce(1)+x2*pce(2);

% Discretization
x0 = -5500;
xend = 1000;%did 1000
dx = 50;
x = x0:dx:xend;
T = 1200; %1200

% Initial density
dis = [0.8 0.2]./pce;
pt = [-2000 0.5*kc;-2000 kjam;0 kjam;0 0];
ux0 = Up(x,pt)'*dis;

%% Calculate density
[UCUP,U1CUP,U2CUP,tCUP] = CU3(x,T,ux0,v,dv,q,xT);
[ULFP,U1LFP,U2LFP,tLFP] = NLLF3(x,T,ux0,v,dv,q,pce);

% graph
figure()
for i = 1:(T/100)
    hold on 
    plot(x,ULFP(:,ceil(i*length(tLFP)/(T/100))))
    legend();
    xlabel("Distance")
    ylabel("Density")
    title("ConstPCE in LF")
end
figure()
contour(tLFP,x,ULFP,'ShowText','on')

figure()
for i = 1:(T/100)
    hold on 
    plot(x,UCUP(:,ceil(i*length(tCUP)/(T/100))))
    legend();
    xlabel("Distance")
    ylabel("Density")
    title("ConstPCE in CU")
end
figure()
contour(tCUP,x,UCUP,'ShowText','on')

figure()
imagesc(tLFP,x,ULFP)
colorbar()
set(gca, 'XLim', tLFP([1 end]), 'YLim', x([1 end]), 'YDir', 'normal')
title("ConstPCE in LF")
xticks([600 800 1000])

figure()
imagesc(tCUP,x,UCUP)
colorbar()
set(gca, 'XLim', tCUP([1 end]), 'YLim', x([1 end]), 'YDir', 'normal')
title("ConstPCE in CU")
xticks([600 800 1000])

%Capacity drop at x for LF
testPT = 200;
figure();
den = UCUP(testPT,:);
flow = q(U1CUP(testPT,:),UCUP(testPT,:),1)+q(U2CUP(testPT,:),UCUP(testPT,:),2);
plot(den,flow,".");

%Capacity drop at x for CU
testPT = 200;
figure();
den = ULFP(testPT,:);
flow = q(U1LFP(testPT,:),ULFP(testPT,:),1)+q(U2LFP(testPT,:),ULFP(testPT,:),2);
plot(den,flow,".");

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
