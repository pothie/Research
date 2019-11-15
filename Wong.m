% Wong experiment 1
%% 
dx = 3/1000; %5m
vmax = [90;120]; %120km/h
T = 1.5;% hour
D = 2; %2km
x = 0:dx:D; %space grid
k0 = 50; % 50 veh/km
pce = [1 2];

% fundamental relation
v = @(xT,n) vmax(n)*exp(-(xT./k0).^2./2); %xT:Total density
dv = @(xT,n) vmax(n)*exp(-(xT./k0).^2./2).*(-1/k0^2*xT);%.*pce(n);
q = @(x,xT,n) x.*v(xT,n); %x:individual density
xT = @(x1,x2) [x1 x2]*[pce(xT,1);pce(xT,2)]';

% Initial density
ux0 = zeros(length(x),2);

% Calculate density
[U,U1,U2,t] = NLLF(x,T,ux0,v,dv,q,xT);
%[U,U1,U2,t] = NLLF4(x,T,ux0,v,dv,q,xT);

%% Plot

av = (q(U1(301,:),U(301,:),1)+q(U2(301,:),U(301,:),2))./U(301,:);
plot(U(301,:),av,'.')
xlabel('Total Density')
ylabel('Average speed')

figure()
qT = (q(U1(301,:),U(301,:),1)+q(U2(301,:),U(301,:),2));
plot(U(301,:),qT,'.')
xlabel('Density')
ylabel('Flow')

figure()
imagesc(t,x,U)
colorbar()
set(gca, 'XLim', t([1 end]), 'YLim', x([1 end]), 'YDir', 'normal')
xticks([0.5 1 1.125 1.175])
xlabel('time')
ylabel('distance')
title('Wong 2-class density graph')