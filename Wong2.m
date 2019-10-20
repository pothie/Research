% Wong
dt = 1.5e-3/60; % hour
dx = 5/1000; %5m
vmax = [120;90]; %3/100m/s (120km/h)
T = 0.02;%1.5h
D = 2; %2km
x = 0:dx:D; %space grid

% time grid
if floor(T/dt)*dt~=T
    t=[0:dt:T T];
else
    t = 0:dt:T;
end
k0 = 50; % 50 veh/km

% fundamental relation
q = @(x,n) vmax(n)*x.*exp(-(x./k0).^2./2);
v = @(x,n) vmax(n)*exp(-(x./k0).^2./2);
% dv = @(x,)
% dq = @(x) v(x).*(1-x.^2/k0^2);

% Initial density
pt = [0 0;0.1 40;0.9 40;1 0];
ux0 = Up(x,pt);
ut1 = zeros(size(t));
uend1 = zeros(size(t));

U1 = NLLF(x,t,ux0,q,dq,ut1,uend1);
U2 = NLLF(x,t,ux0,q,dq,ut1,uend1);

for i = 1:5
hold on
plot(x,U1(:,300*i));
end
xlabel('distance')
ylabel('density')

legend()
hold off
figure()
for i=1:6
    plot(t,q(U1(40*i,:)),'.')
    hold on
end
xlabel('time');
ylabel('flow');
hold off
[X,Y] = meshgrid(x,t);
contour(Y,X,U1',200)
colorbar()