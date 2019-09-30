% Wong
dt = 1.5e-3/60; %hour seconds
dx = 5/1000; %5m
vmax = [120;90]; %120km/h
T = 1.5;% hour
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
dq = @(x) v(x).*(1-x.^2/k0^2);

% Initial density
u0 = zeros(size(x));
pt = [0.25 10;0.5 50;1 50;1.25 10];
u1 = Up(t,pt);
uend1 = ones(size(t))*0;

% Calculate density
U = NLLF(x,t,u0,q,dq,u1,uend1);
plot(U(:,:),v(U(:,:)),'.');
xlabel('distance')
ylabel('speed')

figure()
plot(U(:,:),q(U(:,:)),'.');
xlabel('density');
ylabel('flow');

figure()
[X,Y] = meshgrid(x,t);
contour(Y,X,U',200)
colorbar()
xlabel('time')
ylabel('distance')
