% Smooth function testing
% Burger's equation
clear all
a=0;
b=2*pi;
q = @(x) x.^2./2; 
dq = @(x) x;
%ux0 = @(x) 0.5+sin(x);
%v = @(x,n) x/2; %xT:Total density
%dv = @(xT,n) 1/2;
%q = @(x,xT,n)  xT./2.*x; %x:individual density
%x0 = 0.2;
T = 0.5;
u_exact = @(x,t) sin(x).*cos(t);
error = zeros(6,1);

for i = 1:4
    m = 10*2^(1+i);
    x = linspace(a,b,m+1)';
    ux0 = sin(x);
    [U,t] = CU(x,T,ux0,q,dq);
    hold on
    plot(x,U(:,end))
    xlabel('x');
    ylabel('u');
    error(i) = norm(u_exact(x,T) - U(:,end),inf);
end

plot(x,u_exact(x,T),'b')
disp(error)
hold off

figure()
loglog(error,'r')
hold on; 
loglog((1/4).^(1:6));

hold off
%
%figure()
%for j = 1:3
%hold on
%plot(x,U(:,j*100));
%legend();
%end

close all
