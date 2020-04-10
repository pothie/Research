%Constant PCE
%clear all
vmax = [35 27.5];%30 27.5
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
dx = 25;
x = x0:dx:xend;
T = 1200; %1200

% Initial density
pt = [-2000 0.5*kc;-2000 kjam;0 kjam;0 0];
d1 = 0.8;
d2 = 1-d1;
d = d1/d2;
dis = [d/(d+pce(2)) 1/(d+pce(2))];
ux0 = Up(x,pt)'*dis;

%% capacity drop for different dis
for i = 1:3
    d1 = 0.3+i*0.2;
    d2 = 1-d1;
    d = d1/d2;
    dis = [d/(d+pce(2)) 1/(d+pce(2))];
    ux0 = Up(x,pt)'*dis;
    [UCUP,U1CUP,U2CUP,~] = CU3(x,T,ux0,v,dv,q,xT);
    testPT = 130;
    den = UCUP(testPT,:);
    flow = q(U1CUP(testPT,:),UCUP(testPT,:),1)+q(U2CUP(testPT,:),UCUP(testPT,:),2);
    hold on
    plot(den,flow,".");
    ylabel("Total Flow (veh/s)")
    xlabel("Total Density (veh/m)")
    ylim([0 0.7])
end
legend("Class 1 50%, Class 2 50%","Class 1 70%, Class 2 30%","Class 1 90%, Class 2 10%")
print -depsc Fcdp
%% Calculate density
[UCUP,U1CUP,U2CUP,tCUP] = CU3(x,T,ux0,v,dv,q,xT);
[ULFP,U1LFP,U2LFP,tLFP] = NLLF3(x,T,ux0,v,dv,q,pce);

[UCUN,U1CUN,U2CUN,tCUN] = CU3(x,T,ux0,v,dv,q,xT);
[ULFN,U1LFN,U2LFN,tLFN] = NLLF3(x,T,ux0,v,dv,q,pce);

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
for i = 1:(T/100)
    hold on 
    plot(x,UCUP(:,ceil(i*length(tCUP)/(T/100))))
    legend();
    xlabel("Distance")
    ylabel("Density")
    title("ConstPCE in CU")
end

figure()
imagesc(tLFP,x,ULFP)
colorbar()
set(gca, 'XLim', tLFP([1 end]), 'YLim', x([1 end]), 'YDir', 'normal')
xticks([600 800 1000])

figure()
imagesc(tCUP,x,UCUP)
colorbar()
set(gca, 'XLim', tCUP([1 end]), 'YLim', x([1 end]), 'YDir', 'normal')
xticks([600 800 1000])
xlabel("Time(s)")
ylabel("Distance(m)")

%Capacity drop at x for CU
testPT = 130;
figure();
den = UCUP(testPT,:);
flow = q(U1CUP(testPT,:),UCUP(testPT,:),1)+q(U2CUP(testPT,:),UCUP(testPT,:),2);
plot(den,flow,".");
figure()
sp = flow./den;
plot(den,sp,".")

%Capacity drop at x for LF
testPT = 130;
figure();
den = ULFP(testPT,:);
flow = q(U1LFP(testPT,:),ULFP(testPT,:),1)+q(U2LFP(testPT,:),ULFP(testPT,:),2);
plot(den,flow,".");
ylabel("Total Flow (veh/s)")
xlabel("Total Density (veh/m)")
%print -depsc Fcd
figure()
sp = flow./den;
plot(den,sp,".")
ylabel("Overall Speed (m/s)")
xlabel("Total Density (veh/m)")
%print -depsc Fh

plot(xCU,UCU(:,end),'k*')
hold on
plot(xLF,ULF(:,end),'y*')
legend("CUP","LFP","CU","LF")
