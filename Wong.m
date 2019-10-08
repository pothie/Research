% Wong
%% 
dx = 5/1000; %5m
vmax = [90;120]; %120km/h
T = 1.5;% hour
D = 2; %2km
x = 0:dx:D; %space grid
k0 = 50; % 50 veh/km

% fundamental relation
v = @(xT,n) vmax(n)*exp(-(xT./k0).^2./2); %xT:Total density
dv = @(xT,n) vmax(n)*exp(-(xT./k0).^2./2).*(-1/k0^2*xT);
q = @(x,xT,n) x.*v(xT,n); %x:individual density

% Initial density
ux0 = zeros(length(x),2);

% Calculate density
[U,U1,U2,t] = NLLF2(x,T,ux0,v,dv,q);

%% Plot
plot(U1(301,:),v(U(301,:),1),'b.');
hold on 
plot(U2(301,:),v(U(301,:),2),'k.');
xlabel('density')
ylabel('speed')

figure()
plot(U(301,:),q(U1(301,:),U(301,:),1),'b.');
hold on 
plot(U(301,:),q(U2(301,:),U(301,:),2),'k.');
xlabel('density')
ylabel('flow')

figure()
imagesc(t,x,U)
colorbar()
set(gca, 'XLim', t([1 end]), 'YLim', x([1 end]), 'YDir', 'normal')
xticks([0.5 1 1.125 1.175])
xlabel('time')
ylabel('distance')
title('Wong 2-class density graph')