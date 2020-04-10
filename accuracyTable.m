clear all
vmax = [120;60]; %120km/h
T = 0.01;% hour
xmin = -pi; %2km
xmax = pi;
k0 = 50; % 50 veh/km
pce = [1 1];

% fundamental relation
v = @(xT,n) vmax(n)*exp(-(xT./k0).^2./2); %xT:Total density
dv = @(xT,m,n) vmax(m)*exp(-(xT./k0).^2./2).*(-1/k0^2*xT).*pce(n);
q = @(x,xT,n) x.*v(xT,n); %x:individual density
xT = @(x1,x2) x1*pce(1)+x2*pce(2);

% Initial density
%pt = [0.25 10;0.5 50;1 50;1.25 10];
%pt = [0 0;0.1 40;0.9 40;1 0];
dis = [0.5 0.5]./pce;

% Error
m = 5;%highest power
error = zeros(3,6);
x = linspace(xmin,xmax,10*2^m+1)';
ux0 = 2*sin(x)*dis;%Up(x,pt)'*dis;
[LF_standard,U1LF_standard,U2LF_standard,~] = NLLF3(x,T,ux0,v,dv,q,pce);
[CU_standard,U1CU_standard,U2CU_standard,~] = CU3(x,T,ux0,v,dv,q,xT);

for i = 1:m-1
    x = linspace(xmin,xmax,10*2^i+1)';
    ux0 = 2*sin(x)*dis;%Up(x,pt)'*dis;
    [ULF,U1LF,U2LF,~] = NLLF3(x,T,ux0,v,dv,q,pce);
    [UCU,U1CU,U2CU,~] = CU3(x,T,ux0,v,dv,q,xT);
%     
%     LF = LF_standard(1:2^(m-i):end,1:2^(m-i):end);
%     LF = LF(:,1:min(length(LF(1,:)),length(ULF(1,:))));
%     CU = CU_standard(1:2^(m-i):end,1:2^(m-i):end);
%      CU = CU(:,1:min(length(CU(1,:)),length(UCU(1,:))));
% %     error(i,1) = norm(ULF(:,1:length(LF(1,:)))-LF,1)*1/(10*2^i);
% %     error(i,2) = norm(UCU(:,1:length(CU(1,:)))-CU,1)*1/(10*2^i);
    error(i,1) = norm(ULF(:,end)-LF_standard(1:2^(m-i):end,end),1)*1/(10*2^i);
    error(i,2) = norm(UCU(:,end)-CU_standard(1:2^(m-i):end,end),1)*1/(10*2^i);
    error(i,3) = norm(U1LF(:,end)-U1LF_standard(1:2^(m-i):end,end),1)*1/(10*2^i);
    error(i,4) = norm(U1CU(:,end)-U1CU_standard(1:2^(m-i):end,end),1)*1/(10*2^i);
    error(i,5) = norm(U2LF(:,end)-U2LF_standard(1:2^(m-i):end,end),1)*1/(10*2^i);
    error(i,6) = norm(U2CU(:,end)-U2CU_standard(1:2^(m-i):end,end),1)*1/(10*2^i);
end

error
conrate = log2(error(1:end-1,:)./error(2:end,:))
%%
% plot(x,LF(:,end))
% xlabel("Distance")
% ylabel("Total density")
% %print -depsc LF
% %plot(x,CU(:,end))
% xlabel("Distance")
% ylabel("Total density")
%print -depsc CU
% accuracyTable
% m=5
% error =
%       0.66713      0.17631
%       0.31504     0.052076
%       0.14129     0.020343
%       0.04817    0.0051114
% conrate =
%        1.0824       1.7594
%        1.1569       1.3561
%        1.5524       1.9928
% m=7
% error =
%        0.7034      0.17619
%       0.35032     0.050706
%        0.1774     0.021054
%      0.084695     0.006152
%      0.036371    0.0016872
%      0.012137   0.00040093
% conrate =
%        1.0057       1.7969
%       0.98164       1.2681
%        1.0667       1.7749
%        1.2195       1.8664
%        1.5834       2.0732
%        