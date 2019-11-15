%Dynamic PCE
vmax = [30 27.5];
vc = 25;
L = [6 18]; %L1 was 6, L2 was 18
kjam = 1/L(1);
kc = kjam/6;
w = (vc*kc)/(kjam-kc);
TH = [1 1.5];

% fundamental relation
% uv and udv are universial variables
v = @(xT,n) uv(xT,n,vmax,vc,L(1));
dxT = @(x1,x2,xT,n) udxT(x1,x2,xT,n,L,TH,vmax,vc,kjam,kc);
dv = @(x1,x2,xT,m,n) udv(xT,m,vmax,vc,L(1)).*dxT(x1,x2,xT,n);
q = @(x,xT,n) x.*v(xT,n);
xT = @(x1,x2) uxT(x1,x2,L,TH,vmax,vc,kjam,kc);
pce = @(xT,n) (L(n)+TH(n).*v(xT,n))./(L(1)+TH(1).*v(xT,1));

% Discretization
x0 = -5500;
xend = 1000;
dx = 2.5;
x = x0:dx:xend;
T = 1200; %1200

% Initial density
pt = [-2000 0.5*kc;-2000 kjam;0 kjam;0 0];
ux0 = Up(x,pt)'*[0.8 0.2]./[pce(Up(x,pt),1);pce(Up(x,pt),2)]';

%% Calculate density
[U,U1,U2,t] = NLLF4(x,T,ux0,v,dv,q,xT);
% graph
figure()
imagesc(t,x,U)
colorbar()
set(gca, 'XLim', t([1 end]), 'YLim', x([1 end]), 'YDir', 'normal')
xticks([600 800 1200])

testpt = floor(length(x)/2);
figure()
av = (q(U1(testpt,:),U(testpt,:),1)+q(U2(testpt,:),U(testpt,:),2))./U(testpt,:);
plot(U(testpt,:),av,'.')
xlabel('Total Density')
ylabel('Average speed')

figure()
qT = (q(U1(testpt,:),U(testpt,:),1)+q(U2(testpt,:),U(testpt,:),2));
plot(U(testpt,:),qT,'.')
xlabel('Density')
ylabel('Flow')

figure()
for i = 1:7
    hold on
    plot(x,U(:,i*2e3));
    legend();
end

%% 
l = zeros(250,4);
for i = 1:7096
xTotal = [U1(:,i) U2(:,i)];
pf = (af(1)-xTotal*bf'-((af(1)-xTotal*bf').^2+4*bf(1)*xTotal*af').^0.5)./(-2*bf(1));
pc = (ac(1)-xTotal*bc'-((ac(1)-xTotal*bc').^2+4*bc(1)*xTotal*ac').^0.5)./(-2*bc(1));
pf(abs(imag(pf))>1e-13) = -2; %drop complex solutions
if length(pc(kc<=pc&pc<=kjam+1e-14))+length(pf(0<=pf&pf<kc-1e-17)) ~= length(x1)
    count = count+1; 
    l(count,1)=i;
    l(count,2) = (length(pc(kc<=pc&pc<=kjam+1e-14)));
    l(count,3)=(length(pf(0<=pf&pf<kc-1e-17)));
    l(count,4)=(length(pc(kc<=pc&pc<=kjam+1e-14))+length(pf(0<=pf&pf<kc-1e-17))-length(x1));
end
end

%%
 for i = 1:7095
xTotal = [U1(:,i) U2(:,i)];
pf = (af(1)-xTotal*bf'-((af(1)-xTotal*bf').^2+4*bf(1)*xTotal*af').^0.5)./(-2*bf(1));
pc = (ac(1)-xTotal*bc'-((ac(1)-xTotal*bc').^2+4*bc(1)*xTotal*ac').^0.5)./(-2*bc(1));
pf(abs(imag(pf))>1e-13) = -2; %drop complex solutions
xTnew(kc-1e-9<=pc&pc<=kjam+1e-14,i) = pc(kc-1e-9<=pc&pc<=kjam+1e-14);
xTnew(0<=pf&pf<kc-1e-9,i) = pf(0<=pf&pf<kc-1e-9);if  any(xTnew(:,i)>kjam+1e-14) | any(xTnew(:,i)<0-1e-14)
disp(i);disp(find(any(xTnew(:,i)>kjam+1e-14) | any(xTnew(:,i)<0-1e-14)));
end
end