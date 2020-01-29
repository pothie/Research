%Driverscriptforsolvingthe1DBurgersequationsusinganDGscheme
clear all
%Orderofmethod(m),numberofelements(N)
m=1;N=80;
%Setproblemparameters
xmin=-5500;xmax=1000;
FinalTime=600;CFL=0.1;
vmax = [30 27.5];
vc = 25;
L = [6 18]; %L1 was 6, L2 was 18
kjam = 1/L(1);
kc = kjam/6;
w = (vc*kc)/(kjam-kc);
%Generatemesh
VX=(xmax-xmin)*(0:N)/N+xmin;r=LegendreGL(m);
x=(ones(m+1,1)*VX(1:N)+0.5*(r+1)*(VX(2:N+1)-VX(1:N)))';%x N*3 vector
h=(xmax-xmin)/N;
%Define initial conditions
pce = [1 2];
pt = [-2000 0.5*kc;-2000 kjam;0 kjam;0 0];
dis = [0.8 0.2]./pce;
u1 = dis(1)*Up(x,pt);
u2 = dis(2)*Up(x,pt);

%SolveProblem
vmax = [90;120]; %120km/h
k0 = 50; % 50 veh/km
% fundamental relation
v = @(xT,n) uv(xT,n,vmax,vc,L(1));
dv = @(xT,m,n) udv(xT,m,vmax,vc,L(1))*pce(n);
q = @(x,xT,n) x.*v(xT,n); %x:individual density
xT = @(x1,x2) pce(1)*x1+pce(2)*x2;

[u1,u2]=DG1DSys(x,u1,u2,h,m,N,CFL,FinalTime,q,v,dv,xT);r=r(:);
% x = x';
% xT = @(x) pce*x;
% [u]=DG1DSysTEST(x,[u1 u2]',h,m,N,CFL,FinalTime,q,v,dv,xT);r=r(:);