%Driverscriptforsolvingthe1DBurgersequationsusinganDGscheme
clear all
%Orderofmethod(m),numberofelements(N)
m=1;N=40;
%Setproblemparameters
xmin=0.0;xmax=1.0;
FinalTime=0.5;CFL=0.1;
%Generatemesh
VX=(xmax-xmin)*(0:N)/N+xmin;r=LegendreGL(m);
x=ones(m+1,1)*VX(1:N)+0.5*(r+1)*(VX(2:N+1)-VX(1:N));
h=(xmax-xmin)/N;
%Defineinitialconditions
u=sin(2*pi*x)+0.8;%periodicBCneeded
%u=(1-sign(x-0.2))/2+1;%ConstantBCneeded
%SolveProblem
[u]=BurgersDG1D(x,u,h,m,N,CFL,FinalTime);r=r(:);
