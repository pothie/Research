%Fastlane quadratic equation testing
vmax = [30 27.5];
vc = 25;
L = [5 30]; %L1 was 6, L2 was 18
kjam = 1/L(1);
kc = kjam/6;
T = [1 2]; %T(2) was 1.5
w = (vc*kc)/(kjam-kc);
if w > L(1)/T(1)
    disp('Attention! w>L(1)/T(1)')
end

af = L+T.*vmax;
ac = T*w*kjam;
bf = -T.*((vmax-vc)/kc);
bc = L-T*w;

% Check the exstence of solution
sbf = @(x) bf*x; % x is a column vector
saf = @(x) af*x;
f = @(x,y) (sbf([x(:) y(:)]')-af(1)).^2+4*bf(1)*saf([x(:) y(:)]');
    %Existence of the solution of the quadratic equation of p 
add = @(x,y) x+y; %expecting k1+k2<pc in free flow
x = 0:0.01:0.2;
[X,Y] = meshgrid(x);
contour(X,Y,reshape(f(X,Y),length(x),length(x)),200)
xlabel('k1')
ylabel('k2')
% colorbar
hold on
contour(X,Y,add(X,Y),[kc,kc],'ShowText','on')
hold off
title('b^2-4ac')

% Check solutions
soln = @(x,y) (af(1)-sbf([x(:) y(:)]')-sqrt(f(x,y)))/(-2*bf(1));
solp = @(x,y) (af(1)-sbf([x(:) y(:)]')+sqrt(f(x,y)))/(-2*bf(1));
figure()
x1 = 0:0.001:0.02;
[X1,Y1] = meshgrid(x1);
contour(X1,Y1,reshape(solp(X1,Y1),length(x1),length(x1)),'ShowText','on')% Too large >pc
xlabel('k1')
ylabel('k2')
title('p = -b+sqrt')

figure()
contour(X1,Y1,reshape(soln(X1,Y1),length(x1),length(x1)),'ShowText','on')% Ok
xlabel('k1')
ylabel('k2')
title('p = -b-sqrt')
