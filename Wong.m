% Wong experiment 1
%% 
clear all
dx = 5/1000; %5m
vmax = [120;60]; %120km/h
T = 0.03;% hour
D = 2; %2km
x = 0:dx:D; %space grid
k0 = 50; % 50 veh/km
pce = [1 1];

% fundamental relation
v = @(xT,n) vmax(n)*exp(-(xT./k0).^2./2); %xT:Total density
dv = @(xT,m,n) vmax(m)*exp(-(xT./k0).^2./2).*(-1/k0^2*xT).*pce(n);
q = @(x,xT,n) x.*v(xT,n); %x:individual density
xT = @(x1,x2) [x1 x2]*pce';

% Initial density
%pt = [0.25 10;0.5 50;1 50;1.25 10];
pt = [0 0;0.1 40;0.9 40;1 0];
dis = [0.5 0.5]./pce;
ux0 = Up(x,pt)'*dis;
%ux0 = zeros(length(x),2);

% Calculate density
[ULF,U1LF,U2LF,tLF] = NLLF3(x,T,ux0,v,dv,q,pce);
%[U,U1,U2,t] = NLLF4(x,T,ux0,v,dv,q,xT);
[UCU,U1CU,U2CU,tCU] = CU3(x,T,ux0,v,dv,q,xT);

%% Plot

figure()
hold on
for i = 1:20
qT = q(U1CU(i*20,:),UCU(i*20,:),1)+q(U2CU(i*20,:),UCU(i*20,:),2);
plot(tCU,qT,'.')
xlabel('time')
ylabel('Flow')
title("flow variation with time")
legend();
end

figure()
imagesc(tLF,x,ULF)
colorbar()
set(gca, 'XLim', tLF([1 end]), 'YLim', x([1 end]), 'YDir', 'normal')
xlabel('time')
ylabel('distance')
title('Wong 2-class density graph')

figure()
imagesc(tCU,x,UCU)
colorbar()
set(gca, 'XLim', tCU([1 end]), 'YLim', x([1 end]), 'YDir', 'normal')
xlabel('time')
ylabel('distance')
title('CU:Wong 2-class density graph')

piece = 5;
figure()
for i = 1:piece
    hold on 
    plot(x,ULF(:,ceil(i*length(tLF)/piece)))
    legend();
    xlabel("Distance")
    ylabel("Density")
    title("Wong in LF")
end
figure()
contour(tLF,x,ULF,'ShowText','on')

figure()
for i = 1:piece
    hold on 
    plot(x,UCU(:,ceil(i*length(tCU)/piece)))
    legend();
    xlabel("Distance")
    ylabel("Density")
    title("Wong in CU")
end
figure()
contour(tCU,x,UCU,'ShowText','on')