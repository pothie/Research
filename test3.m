% nonlinear burger's
figure()

T = 1/5;
a = -1;
b = 1;
for i=1:4
    dx=0.001;
    dt=dx*0.1;%LF: dx^2 LW:0.5*dx
    m=1/dx;
    x=linspace(a,b,m+1);
    x0 = 0.2;
    u0=U0(x,x0,2,-1);
    exact = @(x,t) U0(x-1*t,x0,2,-1);
    f = @(x) x.^2;

    if floor(T/dt)*dt~=T
        t=[0:dt:T T];
    else
        t = 0:dt:T;
    end

    u = NLLW(x,t,u0,f);

    error=norm(exact(x,T)-u(:,end),inf);
    subplot(2,2,i);
    plot(x,exact(x,T),'r')
    hold on
    %plot(x,u0,'k')
    plot(x,u(:,end))
    legend('exact','approximation')
    %figure()
    %[X,Y] = meshgrid(x,0:dt:T);
    %contour(X,Y,u',5,'ShowText','on')

end