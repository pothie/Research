%Driverscriptforsolvingthe1DBurgersequationsusinganDGscheme
%clear all
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
%Generatemesh
VX=(xmax-xmin)*(0:N)/N+xmin;r=LegendreGL(m);
x=(ones(m+1,1)*VX(1:N)+0.5*(r+1)*(VX(2:N+1)-VX(1:N)));
% x(:,1) = x(:,1)+1;
% x(:,end) = x(:,end)-1;
h=(xmax-xmin)/N;
%Define initial conditions
pce = [1 L(2)/L(1)];
pt = [-2000 0.5*kc;-2000 kjam;0 kjam;0 0];
d1 = 0.8;
d2 = 1-d1;
d = d1/d2;
dis1 = d./(d+pce(2));
dis2 = 1./(d+pce(2));
u1p = dis1.*Up(x,pt);
u2p = dis2.*Up(x,pt);

%SolveProblem
% fundamental relation
v = @(xT,n) uv(xT,n,vmax,vc,L(1));
dv = @(xT,m,n) udv(xT,m,vmax,vc,L(1))*pce(n);
q = @(x,xT,n) x.*v(xT,n); %x:individual density
xT = @(x1,x2) pce(1)*x1+pce(2)*x2;

% x = x';
% xT = @(x) pce*x;
% [u]=DG1DSysTEST(x,[u1 u2]',h,m,N,CFL,FinalTime,q,v,dv,xT);r=r(:);
[u1p,u2p,tgrid]=DG1DSys(x,u1p,u2p,h,m,N,CFL,FinalTime,q,v,dv,xT);r=r(:);

% Graph
UDGP = xT(u1p,u2p);
%%
plot(x,UDGP,'Color',[0, 0.4470, 0.7410])
ylim([0 0.18])
xlabel("Distance (m)")
ylabel("Total density (veh/m)")

figure()
plot(x,ULFP(:,end))
hold on
plot(x,UCUP(:,end))
VX=(xmax-xmin)*(0:N)/N+xmin;r=LegendreGL(m);
xV=ones(m+1,1)*VX(1:N)+0.5*(r+1)*(VX(2:N+1)-VX(1:N));
plot(xV,UDGP,'k')
xlabel("Distance (m)")
ylabel("Total density (veh/m)")
legend('LF','CU','DG')
