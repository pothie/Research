% Wong
dt = 1.5e-3/60; %seconds
dx = 5/1000; %5m
v1max = 90; %3/100m/s (120km/h)
T = 0.015;%1.5*60*60; % 1.5h
D = 2; %2km
x = 0:dx:D; %space grid
% time grid
if floor(T/dt)*dt~=T
    t=[0:dt:T T];
else
    t = 0:dt:T;
end
k0 = 50; % 50 veh/km
q = @(x) v1max*x.*exp(-(x./k0).^2./2);
v = @(x) v1max*exp(-(x./k0).^2./2);

U = NLLF(x,t,Up(x),q);

plot(x,U(:,end));
xlabel('distance')
ylabel('density')
figure()
for i=1:10
    plot(t,q(U(40*i,:)),'.')
    hold on
end
xlabel('time');
ylabel('flow');
hold off
figure()
%plot(x,q(U(:,end)))
[X,Y] = meshgrid(x,t);
contour(Y,X,U',200)
colorbar()
