%Driverscriptforsolvingthe1DBurgersequationsusinganDGscheme
clear all
%Orderofmethod(m),numberofelements(N)
m=2;N=400;
%Setproblemparameters
xmin=-5500;xmax=1000;
FinalTime=600;CFL=0.1;
vmax = [35 27.5];
vc = 25;
L = [6 18]; %L1 was 6, L2 was 18
kjam = 1/L(1);
kc = kjam/6;
w = (vc*kc)/(kjam-kc);
TH = [1 1.5];
%Generatemesh
VX=(xmax-xmin)*(0:N)/N+xmin;r=LegendreGL(m);
x=(ones(m+1,1)*VX(1:N)+0.5*(r+1)*(VX(2:N+1)-VX(1:N)));
h=(xmax-xmin)/N;
%SolveProblem
% fundamental relation
v = @(xT,n) uv(xT,n,vmax,vc,L(1));
dxT = @(x1,x2,xT,n) udxT(x1',x2',xT,n,L,TH,vmax,vc,kjam,kc);
dv = @(x1,x2,xT,m,n) udv(xT,m,vmax,vc,L(1)).*dxT(x1',x2',xT,n);
q = @(x,xT,n) x.*v(xT,n);
xT = @(x1,x2) uxT(x1,x2,L,TH,vmax,vc,kjam,kc);%split x
pce = @(xT,n) (L(n)+TH(n).*v(xT,n))./(L(1)+TH(1).*v(xT,1));
%Define initial conditions
pt = [-2000 0.5*kc;-2000 kjam;0 kjam;0 0];
d1 = 0.8;
d2 = 1-d1;
d = d1/d2;
dis1 = d./(d+pce(Up(x,pt),2));
dis2 = 1./(d+pce(Up(x,pt),2));
u1 = dis1.*Up(x,pt);
u2 = dis2.*Up(x,pt);

[u1,u2,tgrid]=DG1DSysDy(x,u1,u2,h,m,N,CFL,FinalTime,q,v,dv,xT);r=r(:);

%%
UDG = xT(u1,u2);
figure()
hold on
plot(x(1,:),UDG(1,:),'b');
plot(x(1,:),UDGP(1,:),'r');
plot(x(1,:),UDGN(1,:),'k');
legend('Dynamic PCE','Constant PCE','No PCE')
ylim([0 0.18])
xlabel("Distance (m)")
ylabel("Total density (veh/m)")
print -depsc FD600DG

figure()
plot(x(1,:),u1(1,:),'b');
hold on
plot(x(1,:),u1p(1,:),'r');
plot(x(1,:),u1n(1,:),'k');
legend('Dynamic PCE','Constant PCE','No PCE')
ylim([0 0.18])
xlabel("Distance (m)")
ylabel("Total density (veh/m)")
%print -depsc FD600DG1

figure()
plot(x(1,:),u2(1,:),'b');
hold on
plot(x(1,:),u2p(1,:),'r');
plot(x(1,:),u2n(1,:),'k');
legend('Dynamic PCE','Constant PCE','No PCE')
ylim([0 0.18])
xlabel("Distance (m)")
ylabel("Total density (veh/m)")
%print -depsc FD600DG2
