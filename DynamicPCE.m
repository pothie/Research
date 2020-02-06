%Dynamic PCE
clear all
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
x0 = -5500;
xend = 1000;
dx = 25;
xCU = x0:dx:xend;
xLF = x0:dx:xend;
T = 1200; %1200

% Initial density
pt = [-2000 0.5*kc;-2000 kjam;0 kjam;0 0];
dis = [0.8 0.2];
ux0CU = Up(xCU,pt)'*dis./[pce(Up(xCU,pt),1);pce(Up(xCU,pt),2)]';
ux0LF = Up(xLF,pt)'*dis./[pce(Up(xLF,pt),1);pce(Up(xLF,pt),2)]';

%% Calculate density
%[U,U1,U2,t] = CU3(x,T,ux0,v,dv,q,xT);
[UCU,U1CU,U2CU,tCU] = CU4(xCU,T,ux0CU,v,dv,q,xT);
[ULF,U1LF,U2LF,tLF] = NLLF4(xLF,T,ux0LF,v,dv,q,xT);

% graph
figure()
for i = 1:(T/100)
    hold on 
    plot(xLF,ULF(:,ceil(i*length(tLF)/(T/100))))
    legend();
    xlabel("Distance")
    ylabel("Density")
    title("DynamicPCE in LF")
end
figure()
contour(tLF,xLF,ULF,'ShowText','on')

figure()
for i = 1:(T/100)
    hold on 
    plot(xCU,UCU(:,ceil(i*length(tCU)/(T/100))))
    legend();
    xlabel("Distance")
    ylabel("Density")
    title("DynamicPCE in CU")
end
figure()
contour(tCU,xCU,UCU,'ShowText','on')

figure()
imagesc(tCU,xCU,UCU)
colorbar()
set(gca, 'XLim', tCU([1 end]), 'YLim', xCU([1 end]), 'YDir', 'normal')
xticks([600 800 1000])

figure()
imagesc(tLF,xLF,ULF)
colorbar()
set(gca, 'XLim', tLF([1 end]), 'YLim', xLF([1 end]), 'YDir', 'normal')
xticks([600 800 1000])

figure()
hold on
for i = 1:10
qT = q(U1CU(i*14,:),UCU(i*14,:),1)+q(U2CU(i*14,:),UCU(i*14,:),2);
plot(tCU,qT,'.')
xlabel('time')
ylabel('Flow')
title("flow variation with time")
legend();
end

%Capacity drop at x for LF
testPT = 200;
figure();
den = UCU(testPT,:);
flow = q(U1CU(testPT,:),UCU(testPT,:),1)+q(U2CU(testPT,:),UCU(testPT,:),2);
plot(den,flow,".");

%Capacity drop at x for CU
testPT = 200;
figure();
den = ULF(testPT,:);
flow = q(U1LF(testPT,:),ULF(testPT,:),1)+q(U2LF(testPT,:),ULF(testPT,:),2);
plot(den,flow,".");

figure()
for i = 1:7
    hold on
    plot(xCU,ULF(:,i*1e2));
    legend();
end
