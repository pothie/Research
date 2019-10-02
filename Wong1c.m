% Wong
% Parameters
dx = 5/1000; %5m
vmax = 90; %120km/h
T = 1.5;% hour
D = 2; %2km
x = 0:dx:D; %space grid
k0 = 50; % 50 veh/km

% fundamental relation
q = @(x) vmax*x.*exp(-(x./k0).^2./2);
v = @(x) vmax*exp(-(x./k0).^2./2);
dq = @(x) v(x).*(1-x.^2/k0^2); %by hand

% Initial density when x=0
u0 = zeros(size(x));

% Calculate density
[U,t] = NLLF(x,T,u0,q,dq);

%Graphs

%plot(U(:,:),v(U(:,:)),'b.');
%xlabel('density')
%ylabel('speed')

%figure()
%plot(U(:,:),q(U(:,:)),'b.');
%xlabel('density');
%ylabel('flow');

%figure()
imagesc(t,x,U)
colorbar()
set(gca, 'XLim', t([1 end]), 'YLim', x([1 end]), 'YDir', 'normal')
xlabel('time')
ylabel('distance')
title('Wong 1-class density graph')