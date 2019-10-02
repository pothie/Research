% Wong
% dt = 1.5e-3/60; %hour
dx = 5/1000; %5m
vmax = [90;120]; %120km/h
T = 0.5;% hour
D = 2; %2km
x = 0:dx:D; %space grid
k0 = 50; % 50 veh/km

% fundamental relation
v = @(xT,n) vmax(n)*exp(-(xT./k0).^2./2); %xT:Total density
dv = @(x,n) v(x,n).*(-1/k0^2*x);
q = @(x,xT,n) x.*v(xT,n); %x:individual density

% Initial density
u0 = zeros(length(x),2);
%pt = [0.25 10;0.5 50;1 50;1.25 10];
%u1 = Up(t,pt);
%uend1 = ones(size(t))*0;

% Calculate density
[U,U1,U2,t] = NLLF2(x,T,u0,v,dv,q);
%plot(U(:,:),v(U(:,:),1),'b.');
%hold on 
%plot(U(:,:),v(U(:,:),2),'k.');
%xlabel('density')
%ylabel('speed')

%figure()
%plot(U1(301,:),q(U1(301,:),1),'b.');
%hold on 
%plot(U2(301,:),q(U2(301,:),2),'k.');
%xlabel('density')
%ylabel('flow')

figure()
imagesc(t,x,U)
colorbar()
set(gca, 'XLim', t([1 end]), 'YLim', x([1 end]), 'YDir', 'normal')
xlabel('time')
ylabel('distance')
title('Wong 2-class density graph')