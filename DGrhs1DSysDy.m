function[rhsu1,rhsu2]=DGrhs1DSysDy(x,u1,u2,h,k,m,N,Ma,S,VtoE,maxvel,f,time,xT)
%function[rhsu]=BurgersDGrhs1D(x,u,h,k,m,N,Ma,S,VtoE,maxvel)
%Purpose:Evaluate the RHS of Burgers equations usinga DGmethod
Imat=eye(m+1);%ue=zeros(N+2,2); 
uT = xT(u1,u2);
%calculate flux
[~,n] = size(u1);
flux1 = zeros(1,1+n);
flux2 = zeros(size(flux1));

flux1(2:end-1) = (f(u1(end,1:end-1),uT(end,1:end-1),1)+f(u1(1,2:end),uT(1,2:end),1))/2 ...
        -maxvel/2.*(u1(1,2:end)-u1(end,1:end-1));
flux2(2:end-1) = (f(u2(end,1:end-1),uT(end,1:end-1),2)+f(u2(1,2:end),uT(1,2:end),2))/2 ...
        -maxvel/2.*(u2(1,2:end)-u2(end,1:end-1));

flux1(1) = f(u1(1,1),uT(1,1),1);
flux1(end) = f(u1(end,end),uT(end,end),1);

flux2(1) = f(u2(1,1),uT(1,1),2);
flux2(end) = f(u2(end,end),uT(end,end),2);

ru1=S'*f(u1,uT,1)-(Imat(:,m+1)*flux1(2:end)-Imat(:,1)*flux1(1:end-1));
ru2=S'*f(u2,uT,2)-(Imat(:,m+1)*flux2(2:end)-Imat(:,1)*flux2(1:end-1));
rhsu1=(h/2*Ma)\ru1;% row
rhsu2=(h/2*Ma)\ru2;
return