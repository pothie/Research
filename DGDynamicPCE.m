%Driverscriptforsolvingthe1DBurgersequationsusinganDGscheme
clear all
%Orderofmethod(m),numberofelements(N)
m=2;N=100;
%Setproblemparameters
xmin=-5500;xmax=1000;
FinalTime=1200;CFL=0.1;
vmax = [30 27.5];
vc = 25;
L = [6 6]; %L1 was 6, L2 was 18
kjam = 1/L(1);
kc = kjam/6;
w = (vc*kc)/(kjam-kc);
%pce = [1 L(2)/(L(1))];
TH = [1 1];
%Generatemesh
VX=(xmax-xmin)*(0:N)/N+xmin;r=LegendreGL(m);
x=(ones(m+1,1)*VX(1:N)+0.5*(r+1)*(VX(2:N+1)-VX(1:N)))';
h=(xmax-xmin)/N;
%SolveProblem
% fundamental relation
v = @(xT,n) uv(xT,n,vmax,vc,L(1));
dxT = @(x1,x2,xT,n) udxT(x1',x2',xT,n,L,TH,vmax,vc,kjam,kc);
dv = @(x1,x2,xT,m,n) udv(xT,m,vmax,vc,L(1)).*dxT(x1',x2',xT,n);
q = @(x,xT,n) x.*v(xT,n);
xT = @(x1,x2) uxT(x1,x2,L,TH,vmax,vc,kjam,kc);%split x
pce = @(xT,n) (L(n)+TH(n).*v(xT,n))./(L(1)+TH(1).*v(xT,1));
%Define initial conditions
pt = [-2000 0.5*kc;-2000 kjam;0 kjam;0 0];
dis1 = 0.5*ones(size(x))./pce(Up(x,pt),1);
dis2 = 0.5*ones(size(x))./pce(Up(x,pt),2);
u1 = dis1.*Up(x,pt);
u2 = dis2.*Up(x,pt);

[u1r,u1l,u2r,u2l,u1m,u2m,tgrid]=DG1DSysDy(x',u1,u2,h,m,N,CFL,FinalTime,q,v,dv,xT);r=r(:);

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
    title("Dy in DG, m=1, N="+num2str(N))
end
else 
for i = 1:piece
    hold on 
    plot(x',[UDGr(:,ceil(i*length(tgrid)/5)) UDGm(:,ceil(i*length(tgrid)/5))...
        UDGl(:,ceil(i*length(tgrid)/5))]')
    xlabel("Distance")
    ylabel("Density")
    title("Dy in DG, m=2, N="+num2str(N))
end
end

figure()
imagesc(tgrid,x(:,1),UDGr)
colorbar()
set(gca, 'XLim', tgrid([1 end]), 'YLim', x([1 end],1), 'YDir', 'normal')
xlabel('time')
ylabel('distance')
title('DG: left values')
