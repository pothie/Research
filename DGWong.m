%Driverscriptforsolvingthe1DBurgersequationsusinganDGscheme
clear all
% Wong
%Orderofmethod(m),numberofelements(N)
m=1;N=400;
%Setproblemparameters
xmin=0.0;xmax=2;
FinalTime=0.03;CFL=0.1;
%Generatemesh
VX=(xmax-xmin)*(0:N)/N+xmin;r=LegendreGL(m);
x=ones(m+1,1)*VX(1:N)+0.5*(r+1)*(VX(2:N+1)-VX(1:N));
h=(xmax-xmin)/N;
%Define initial conditions
pce = [1 1];
pt = [0 0;0.1 40;0.9 40;1 0];
dis = [0.5 0.5]./pce;
u1 = dis(1)*Up(x,pt);
u2 = dis(2)*Up(x,pt);
%SolveProblem
vmax = [120;60]; %120km/h
k0 = 50; % 50 veh/km
% fundamental relation
v = @(xT,n) vmax(n)*exp(-(xT./k0).^2./2); %xT:Total density
dv = @(xT,m,n) vmax(m)*exp(-(xT./k0).^2./2).*(-1/k0^2*xT).*pce(n);
q = @(x,xT,n) x.*v(xT,n); %x:individual density

% u = [u1;u2];
% xT = @(x) pce*x;
% [u]=DG1DSysTEST(x,u,h,m,N,CFL,FinalTime,q,v,dv,xT);r=r(:);
xT = @(x1,x2) pce(1)*x1+pce(2)*x2;
[u1r,u1l,u2r,u2l,tgrid]=DG1DSys(x',u1',u2',h,m,N,CFL,FinalTime,q,v,dv,xT);r=r(:);