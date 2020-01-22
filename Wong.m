% Wong experiment 1
%% 
dx = 5/1000; %5m
vmax = [90;120]; %120km/h
T = 0.1;% hour
D = 2; %2km
x = 0:dx:D; %space grid
k0 = 50; % 50 veh/km
pce = [1 2];

% fundamental relation
v = @(xT,n) vmax(n)*exp(-(xT./k0).^2./2); %xT:Total density
dv = @(xT,n,m) vmax(n)*exp(-(xT./k0).^2./2).*(-1/k0^2*xT).*pce(n);
q = @(x,xT,n) x.*v(xT,n); %x:individual density
xT = @(x1,x2) [x1 x2]*pce';

% Initial density
pt = [0.25 10;0.5 50;1 50;1.25 10];
dis = [0.8 0.2]./pce;
ux0 = Up(x,pt)'*dis;
%ux0 = zeros(length(x),2);

% Calculate density
[ULF,U1LF,U2LF,tLF] = NLLF3(x,T,ux0,v,dv,q,pce);
%[U,U1,U2,t] = NLLF4(x,T,ux0,v,dv,q,xT);
[UCU,U1CU,U2CU,tCU] = CU3(x,T,ux0,v,dv,q,pce);

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