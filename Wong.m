% Wong experiment 1
%% 
clear all
%dx = 5/1000; %5m
vmax = [120;60]; %120km/h
T = 0.01;% hour
D = 2; %2km
x = linspace(0,D,400); %space grid
k0 = 50; % 50 veh/kmcl
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
%figure()
%% Plot

%Capacity drop at testPT for CU
testPT = 100;%ceil(length(x)/2);
figure();
den = UCU(testPT,:);
flow = q(U1CU(testPT,:),UCU(testPT,:),1)+q(U2CU(testPT,:),UCU(testPT,:),2);
plot(den,flow,".");

%Capacity drop at x for LF
testPT = 300;
figure();
den = ULF(testPT,:);
flow = q(U1LF(testPT,:),ULF(testPT,:),1)+q(U2LF(testPT,:),ULF(testPT,:),2);
plot(den,flow,".");
%Hysteresis for LF
figure()
sp = flow./den;
plot(den,sp,".");
%same graph for LF and Cu

figure()
imagesc(tLF,x,ULF)
colorbar()
set(gca, 'XLim', tLF([1 end]), 'YLim', x([1 end]), 'YDir', 'normal')
xlabel('Time (h)')
ylabel('Distance (km)')
%title('Wong 2-class density graph')

figure()
imagesc(tCU,x,UCU)
colorbar()
set(gca, 'XLim', tCU([1 end]), 'YLim', x([1 end]), 'YDir', 'normal')
xlabel('Time (h)')
ylabel('Distance (km)')
%title('CU:Wong 2-class density graph')

piece = 3;
figure()
for i = 1:piece
    hold on 
    plot(x,ULF(:,ceil(0.5*i*length(tLF)/piece)))
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
    plot(x,UCU(:,ceil(0.5*i*length(tCU)/piece)))
    xlabel("Distance (km)")
    ylabel("Density (veh/km)")
end
figure()
contour(tCU,x,UCU,'ShowText','on')

figure()
plot(x,UCU(:,ceil(length(tCU)/4)))
ylim([0 45])
xlabel("Distance (km)")
ylabel("Total density (veh/km)")

figure()
plot(x,ULF(:,ceil(length(tLF)/4)))
ylim([0 45])
xlabel("Distance (km)")
ylabel("Total density (veh/km)")