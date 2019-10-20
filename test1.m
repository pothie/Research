% Smooth function testing
% Burger's equation
a=0;
b=2*pi;
%q = @(x) x.^2./2; 
%dq = @(x) x;
v = @(x,n) x/2; %xT:Total density
dv = @(xT,n) 1/2;
q = @(x,xT,n)  xT./2.*x; %x:individual density
%x0 = 0.2;
T = 0.5;
u_exact = @(x,t) sin(x).*cos(t);
error = zeros(6,1);

for i =1:6
    m = 2^(5+i);
    x = linspace(a,b,m+1);
    ux0 = sin(x)';
    [U,U1,U2,t] = NLLF2test(x,T,ux0,v,dv,q);
    
    hold on
    plot(x,U(:,end))
    xlabel('x');
    ylabel('u');
    error(i) = norm(u_exact(x,T)-U(:,end)');
end
plot(x,u_exact(x,T),'b')
disp(error)
hold off

figure()
loglog(error,'r')
hold on; 
loglog((1/2).^(1:6));

