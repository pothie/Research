% Discontinuous function testing
% Burger's equation
a=0;
b=2*pi;
q = @(x) x.^2./2; 
dq = @(x) x;
% kjam = 1/6;
% kc = kjam/6;
% pt = [-21 0.5*kc;-20+1 kjam;0-1 kjam;0+1 0];
%ux0 = @(x) 0.5+sin(x);
%v = @(x,n) x/2; %xT:Total density
%dv = @(xT,n) 1/2;
%q = @(x,xT,n)  xT./2.*x; %x:individual density
%x0 = 0.2;
T = 0.5;
u_exact = @(x,t) sin(x).*cos(t);
%error = zeros(6,1);

for i = 1:4
    m(i) = 10*2^(1+i);
    x = linspace(a,b,m(i)+1)';
    ux0 = sin(x);
    [U,t] = NLLF(x,T,ux0,q,dq);
    hold on
    plot(x,U(:,end))
    xlabel('x');
    ylabel('u');
    legend();
    error(i) = norm(u_exact(x,T) - U(:,end),inf)*1/m(i);
end

plot(x,u_exact(x,T),'b')
disp(error)
hold off

figure()
loglog(error,'r')
hold on; 
loglog((1/2).^(1:6));

hold off

%figure()
%for j = 1:3
%hold on
%plot(x,U(:,j*100));
%legend();
%end
T = 20;
m = 160;
x = linspace(a,b,m+1)';
ux0 = Up(x,pt)';
[U1,~] = CU(x,T,ux0,q,dq);
[U2,~] = NLLF(x,T,ux0,q,dq);
ULF = U2(:,end);
UCU = U1(:,end);
plot(x,ULF,'rx');
hold on;
plot(x,UCU,'b.');
plot(x,ux0');
%plot(x,u_exact(x,T),'b')
xlabel('x');
legend('Lax-Friedrichs','Central Upwind');%,'exact solution');
ylabel('u');

