%Driverscriptforsolvingthe1DBurgersequationsusinganDGscheme
clear all
% Wong
%Order of method(m),number of elements(N)
n=5;
m=2;N=400;
%Set problem parameters
xmin=0;xmax=2;
CFL=0.01;
%Generate mesh
VX=(xmax-xmin)*(0:N)/N+xmin;r=LegendreGL(m);
x=ones(m+1,1)*VX(1:N)+0.5*(r+1)*(VX(2:N+1)-VX(1:N));
h=(xmax-xmin)/N;
FinalTime=0.01;
%Define initial conditions
pce = [1 1];
pt = [0 0;0.1 40;0.9 40;1 0];
dis = [0.5 0.5]./pce;
u1 = dis(1)*Up(x,pt);
u2 = dis(2)*Up(x,pt);

%Solve Problem
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

[u1,u2,~]=DG1DSys(x,u1,u2,h,m,N,CFL,FinalTime,q,v,dv,xT);%r=r(:);
UDG = xT(u1,u2);

% for i = 1:n-1
%     N = 10*2^i;
%     VX=(xmax-xmin)*(0:N)/N+xmin;r=LegendreGL(m);
%     x=ones(m+1,1)*VX(1:N)+0.5*(r+1)*(VX(2:N+1)-VX(1:N));
%     h=(xmax-xmin)/N;
%     u1 = 2*dis(1)*sin(x);%*Up(x,pt);
%     u2 = 2*dis(2)*sin(x);%Up(x,pt);
%     [u1,u2,tgrid]=DG1DSys(x,u1,u2,h,m,N,CFL,FinalTime,q,v,dv,xT);%r=r(:);
%     UDG = xT(u1,u2);
%     UDG = UDG([1 end],:);
%     DG = zeros(2,N);
%     DG(1,:) = UDG_standard(1,1:2^(n-i):end);
%     DG(end,:) = UDG_standard(end,2^(n-i):2^(n-i):end);
%    
%     error(i,1) = norm(DG-UDG,1)*1/(10*2^i);
% end
% m
% error
% log2(error(1:end-1,:)./error(2:end,:))

%Graph
piece = 1;
figure()
plot(x,UDG,'Color',[0, 0.4470, 0.7410])
ylim([0 45])
xlabel("Distance (km)")
ylabel("Total density (veh/km)")



