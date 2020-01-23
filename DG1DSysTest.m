%Driverscriptforsolvingthe1DBurgersequationsusinganDGscheme
clear all
%Orderofmethod(m),numberofelements(N)
m=1;N=40;
%Setproblemparameters
xmin=0.0;xmax=2;
FinalTime=0.05;CFL=0.1;
%Generatemesh
VX=(xmax-xmin)*(0:N)/N+xmin;r=LegendreGL(m);
x=ones(m+1,1)*VX(1:N)+0.5*(r+1)*(VX(2:N+1)-VX(1:N));
h=(xmax-xmin)/N;
%Define initial conditions
pce = [1 2];
pt = [0.25 10;0.5 50;1 50;1.25 10];
dis = [0.8 0.2]./pce;
u(1:2,:) = dis(1)*Up(x,pt);
u(3:4,:) = dis(2)*Up(x,pt);
%SolveProblem
vmax = [90;120]; %120km/h
k0 = 50; % 50 veh/km
% fundamental relation
v = @(xT,n) vmax(n)*exp(-(xT./k0).^2./2); %xT:Total density
dv = @(xT,n,m) vmax(n)*exp(-(xT./k0).^2./2).*(-1/k0^2*xT).*pce(n);
q = @(x,xT,n) x.*v(xT,n); %x:individual density
xT = @(x) pce*x;

[u]=DG1DSys(x,u,h,m,N,CFL,FinalTime,q,v,dv,xT);r=r(:);
