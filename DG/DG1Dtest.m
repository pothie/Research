%Driverscriptforsolvingthe1DBurgersequationsusinganDGscheme
clear all
%Orderofmethod(m),numberofelements(N)
n=7;
m=3;N=10*2^n;
%Setproblemparameters
xmin=-pi;xmax=pi;
FinalTime=0.015;CFL=0.1;
%Generatemesh
VX=(xmax-xmin)*(0:N)/N+xmin;r=LegendreGL(m);
x=ones(m+1,1)*VX(1:N)+0.5*(r+1)*(VX(2:N+1)-VX(1:N));
h=(xmax-xmin)/N;
%Define initial conditions
%u=sin(2*pi*x)+0.8;%periodicBCneeded
%pt = [0.25 10;0.5 50;1 50;1.25 10];
u = 10*sin(x);%Up(x,pt);
%u=(1-sign(x-0.2))/2+1;%ConstantBCneeded
%SolveProblem
vmax = 90; %120km/h
k0 = 50; % 50 veh/km 
q = @(x) vmax*x.*exp(-(x./k0).^2./2);
v = @(x) vmax*exp(-(x./k0).^2./2);
dq = @(x) v(x).*(1-x.^2/k0^2); %by hand
% q = @(x) x.^2;
% dq = @(x) 2*x;
[DG_standard]=BurgersDG1D(x,u,h,m,N,CFL,FinalTime,q,dq);r=r(:);

for i=1:n-1
    N = 10*2^i;
    VX=(xmax-xmin)*(0:N)/N+xmin;r=LegendreGL(m);
    x=ones(m+1,1)*VX(1:N)+0.5*(r+1)*(VX(2:N+1)-VX(1:N));
    h=(xmax-xmin)/N;
    u = 10*sin(x);%Up(x,pt);
    [uDG]=BurgersDG1D(x,u,h,m,N,CFL,FinalTime,q,dq);r=r(:);
    DG = DG_standard(:,1:2^(n-i):end);
    DG = DG(:,1:min(length(DG(1,:)),length(uDG(1,:))));
    error(i) = norm(uDG(:,min(length(DG(1,:)),length(uDG(1,:))))-DG(:,end),1)*1/(10*2^i);
    plot(x,uDG)
    hold on
end
m
error'
log2(error(1:end-1)./error(2:end))'
