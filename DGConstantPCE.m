%Driverscriptforsolvingthe1DBurgersequationsusinganDGscheme
clear all
%Orderofmethod(m),numberofelements(N)
m=1;N=250;
%Setproblemparameters
xmin=-5500;xmax=1000;
FinalTime=1200;CFL=0.1;
vmax = [30 27.5];
vc = 25;
L = [6 18]; %L1 was 6, L2 was 18
kjam = 1/L(1);
kc = kjam/6;
w = (vc*kc)/(kjam-kc);
%Generatemesh
VX=(xmax-xmin)*(0:N)/N+xmin;r=LegendreGL(m);
x=(ones(m+1,1)*VX(1:N)+0.5*(r+1)*(VX(2:N+1)-VX(1:N)))';%x N*3 vector
x(:,1) = x(:,1)+1;
x(:,end) = x(:,end)-1;
h=(xmax-xmin)/N;
%Define initial conditions
pce = [1 L(2)/L(1)];
pt = [-2000 0.5*kc;-2000 kjam;0 kjam;0 0];
dis = [0.8 0.2]./pce; %same?
u1 = dis(1)*Up(x,pt);
u2 = dis(2)*Up(x,pt);

%SolveProblem
% fundamental relation
v = @(xT,n) uv(xT,n,vmax,vc,L(1));
dv = @(xT,m,n) udv(xT,m,vmax,vc,L(1))*pce(n);
q = @(x,xT,n) x.*v(xT,n); %x:individual density
xT = @(x1,x2) pce(1)*x1+pce(2)*x2;

% x = x';
% xT = @(x) pce*x;
% [u]=DG1DSysTEST(x,[u1 u2]',h,m,N,CFL,FinalTime,q,v,dv,xT);r=r(:);
[u1r,u1l,u2r,u2l,u1m,u2m,tgrid]=DG1DSys(x',u1,u2,h,m,N,CFL,FinalTime,q,v,dv,xT);r=r(:);

% Graph
UDGr = xT(u1r,u2r);
UDGl = xT(u1l,u2l);
UDGm = xT(u1m,u2m);
piece = 5;
figure()
if m==1
for i = 1:piece
    hold on 
    plot(x',[UDGr(:,ceil(i*length(tgrid)/5)) UDGl(:,ceil(i*length(tgrid)/5))]');
    xlabel("Distance")
    ylabel("Density")
    title("CP in DG, m=1, N="+num2str(N))
end
else 
for i = 1:piece
    hold on 
    plot(x',[UDGr(:,ceil(i*length(tgrid)/5)) UDGm(:,ceil(i*length(tgrid)/5))...
        UDGl(:,ceil(i*length(tgrid)/5))]')
    xlabel("Distance")
    ylabel("Density")
    title("CP in DG, m=2, N="+num2str(N))
end
end

figure()
imagesc(tgrid,x(:,1),UDGl)
colorbar()
set(gca, 'XLim', tgrid([1 end]), 'YLim', x([1 end],1), 'YDir', 'normal')
xlabel('time')
ylabel('distance')
title('DG: left values')
