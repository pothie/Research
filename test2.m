% Non linear testing
% Burger's equation
a=0;
b=1;
q = @(x) x.^2; 
x0 = 0.2;
u_exact = @(x,t) U0(x-3*t,x0,2,1);
T = 1/5;
error = zeros(4,1);

for i =1:6
    m = 2^(5+i);
    x = linspace(a,b,m+1);
    u0 = U0(x,x0,2,1);

    %grid/mesh size
    dx=(b-a)/m;
    dt=dx*0.05;
    n_t=floor(T/dt);

    if n_t*dt~=T
        t=[0:dt:T T];
    else
        t = 0:dt:T;
    end
    u = NLLF(x,t,u0,q);

    %[X,Y] = meshgrid(x,(1:n_t+1)*dt);
    %contour(X,Y,u','ShowText','on')
    %surf(X,Y,u')
    %figure()
    subplot(2,3,i);
    hold on
    plot(x, u(:,end), 'b-', 'LineWidth', 2);
    %plot(x,u0);
    plot(x,u_exact(x,T),'r');
    error(i)=norm(u_exact(x,T)-u(:,end)',Inf);
    xlabel('x');
    ylabel('u');
end
disp(error)
hold off


