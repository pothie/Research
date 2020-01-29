%Driverscriptforsolvingthe1DBurgersequationsusinganDGscheme
% clear all
%Orderofmethod(m),numberofelements(N)
m=1;N=40;
%Setproblemparameters
xmin=-5500;xmax=1000;
FinalTime=600;CFL=0.1;
vmax = [30 27.5];
vc = 25;
L = [6 18]; %L1 was 6, L2 was 18
kjam = 1/L(1);
kc = kjam/6;
w = (vc*kc)/(kjam-kc);
%pce = [1 L(2)/(L(1))];
TH = [1 1.5];
%Generatemesh
VX=(xmax-xmin)*(0:N)/N+xmin;r=LegendreGL(m);
x=ones(m+1,1)*VX(1:N)+0.5*(r+1)*(VX(2:N+1)-VX(1:N));
h=(xmax-xmin)/N;
%SolveProblem
% fundamental relation
v = @(xT,n) uv(xT',n,vmax,vc,L(1));
dxT = @(x1,x2,xT,n) udxT(x1',x2',xT',n,L,TH,vmax,vc,kjam,kc);
dv = @(x1,x2,xT,m,n) udv(xT',m,vmax,vc,L(1)).*dxT(x1',x2',xT',n);
q = @(x,xT,n) x'.*v(xT,n);
xT = @(x) uxT(x(1,:)',x(2,:)',L,TH,vmax,vc,kjam,kc);%split x
pce = @(xT,n) (L(n)+TH(n).*v(xT,n))./(L(1)+TH(1).*v(xT,1));
%Define initial conditions
pt = [-2000 0.5*kc;-2000 kjam;0 kjam;0 0];
dis1 = 0.8*ones(size(x))./pce(Up(x,pt),1)';
dis2 = 0.2*ones(size(x))./pce(Up(x,pt),2)';
u(1:2,:) = dis1.*Up(x,pt);
u(3:4,:) = dis2.*Up(x,pt);

[u]=DG1DSysDy(x,u,h,m,N,CFL,FinalTime,q,v,dv,xT);r=r(:);
