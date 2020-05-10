%Driverscriptforsolvingthe1DwaveequationsusinganDGscheme
clear all
%Orderofmethod(m),numberofelements(N)
n=6;
m=1;
%Setproblemparameters
FinalTime=0.1;
CFL=0.001;
%FinalTime=0.098175*CFL;
xmin=-pi;xmax=pi;
xT = @(x1,x2) x1+x2;

vmax = [120;60];
k0 = 50;
pce = [1 1];
v = @(xT,n) 1;%vmax(n)*exp(-(xT./k0).^2./2); %xT:Total density
dv = @(xT,m,n) 0;%vmax(m)*exp(-(xT./k0).^2./2).*(-1/k0^2*xT).*pce(n);
q = @(x,xT,n) x.*v(xT,n); %x:individual density

for i=1:n
    N = 10*2^i;
    VX=(xmax-xmin)*(0:N)/N+xmin;r=LegendreGL(m);
    x=ones(m+1,1)*VX(1:N)+0.5*(r+1)*(VX(2:N+1)-VX(1:N));
    h=(xmax-xmin)/N;
    u = sin(x);
    u1 = u;
    u2 = u;
    [u1,u2,~]=DG1DSys(x,u1,u2,h,m,N,CFL,FinalTime,q,v,dv,xT);
    [u]=LinwaveDG1D(x,u,h,m,N,CFL,FinalTime);
    error(i,1) = norm(u1-sin(x-FinalTime),inf);%*1/(2^i);
    error(i,2) = norm(u-sin(x-FinalTime),inf);%*1/(2^i);
end
m
error
log2(error(1:end-1,:)./error(2:end,:))

% m=1
% error =
%     0.0084325
%     0.0017264
%    0.00031986
%    5.7381e-05
%     8.323e-06
% conrate =
%        2.2882
%        2.4322
%        2.4788
%        2.7854
