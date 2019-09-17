% Wong
dt = 0.09; %seconds
dx = 5; %m
v1max = 100/3; %m/s (120km/h)
T = 60;%1.5*60*60; % 1.5h
D = 2000; %2km
x = 0:dx:D; %space grid
% time grid
if floor(T/dt)*dt~=T
    t=[0:dt:T T];
else
    t = 0:dt:T;
end
k0 = 50/1000; % 50 veh/km
q = @(x) v1max*x.*exp(-(x./k0).^2./2);

U = NLLF(x,t,U0(x,0,1,0),q);

plot(x,U(:,end));
figure()
plot(t,U(100,:),'.')
figure()
%plot(x,q(U(:,end)))
[X,Y] = meshgrid(x,t);
contour(Y,X,U')
colorbar()




