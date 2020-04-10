%Dynamic PCE
%clear all
vmax = [35 27.5];
vc = 25;
L = [6 18]; %L1 was 6, L2 was 18
kjam = 1/L(1);
kc = kjam/6;
w = (vc*kc)/(kjam-kc);
TH = [1 1.5];%1 1.5

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
d1 = 0.8;
d2 = 1-d1;
d = d1/d2;
disCU = [d./(d+pce(Up(xCU,pt),2));1./(d+pce(Up(xCU,pt),2))]';
disLF = [d./(d+pce(Up(xLF,pt),2));1./(d+pce(Up(xLF,pt),2))]';
ux0CU = repmat(Up(xCU,pt),2,1)'.*disCU;
ux0LF = repmat(Up(xLF,pt),2,1)'.*disLF;

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
for i = 1:(T/100)
    hold on 
    plot(xCU,UCU(:,ceil(i*length(tCU)/(T/100))))
    legend();
    xlabel("Distance")
    ylabel("Density")
    title("DynamicPCE in CU")
end

% capacity drop and hysteresis
testPT = floor(length(xCU)/2);
figure();
den = UCU(testPT,:);
flow = q(U1CU(testPT,:),UCU(testPT,:),1)+q(U2CU(testPT,:),UCU(testPT,:),2);
plot(den,flow,".");
figure()
sp = flow./den;
plot(den,sp,".")

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
for i = 1:7
    hold on
    plot(xCU,ULF(:,i*length(tCU/7)));
    legend();
end


testPT = floor(length(tLF)/2);
plot(x,ULF(:,testPT),'b')
hold on
testPT = floor(length(tLFP)/2);
plot(x,ULFP(:,testPT),'r')
testPT = floor(length(tLFN)/2);
plot(x,ULFN(:,testPT),'k')
ylim([0 0.18])
legend("Dynamic PCE","Constant PCE","No PCE")
ylabel("Total Density (veh/m)")
xlabel("Distance(m)")
print -depsc FD600LF

figure()
testPT = floor(length(tCU)/2);
plot(x,UCU(:,testPT),'b')
hold on
testPT = floor(length(tCUP)/2);
plot(x,UCUP(:,testPT),'r')
testPT = floor(length(tCUN)/2);
plot(x,UCUN(:,testPT),'k')
legend("Dynamic PCE","Constant PCE","No PCE")
ylabel("Total Density (veh/m)")
xlabel("Distance(m)")
print -depsc FD600CU

figure()
testPT = floor(length(tCU)/2);
plot(x,U1CU(:,testPT),'b')
hold on
testPT = floor(length(tCUP)/2);
plot(x,U1CUP(:,testPT),'r')
testPT = floor(length(tCUN)/2);
plot(x,U1CUN(:,testPT),'k')
legend("Dynamic PCE","Constant PCE","No PCE")
ylabel("Density of Class 1(veh/m)")
xlabel("Distance(m)")
print -depsc FD600CU1

figure()
testPT = floor(length(tCU)/2);
plot(x,U2CU(:,testPT),'b')
hold on
testPT = floor(length(tCUP)/2);
plot(x,U2CUP(:,testPT),'r')
testPT = floor(length(tCUN)/2);
plot(x,U2CUN(:,testPT),'k')
legend("Dynamic PCE","Constant PCE","No PCE")
ylabel("Density of Class 2 (veh/m)")
xlabel("Distance(m)")
print -depsc FD600CU2
