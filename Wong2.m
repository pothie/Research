% Wong
dt = 1.5e-3/60; % hour(seconds)
dx = 5/1000; %5m
v1max = 90; %3/100m/s (120km/h)
%for i = 1:5
T = 0.02;%1.5*60*60; % 1.5h
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
q = @(x) v1max*x.*exp(-(x./k0).^2./2);
v = @(x) v1max*exp(-(x./k0).^2./2);
dq = @(x) v(x).*(1-x.^2/k0^2);

% Initial density
pt = [0 0;0.1 40;0.9 40;1 0];
u0 = Up(x,pt);
u1 = zeros(size(t));
uend1 = zeros(size(t));

U = NLLF(x,t,u0,q,dq,u1,uend1);

for i = 1:5
hold on
plot(x,U(:,300*i));
end
xlabel('distance')
ylabel('density')

legend()
hold off
figure()
for i=1:6
    plot(t,q(U(40*i,:)),'.')
    hold on
end
xlabel('time');
ylabel('flow');
hold off
[X,Y] = meshgrid(x,t);
contour(Y,X,U',200)
colorbar()