function[rhsu]=DGrhs1DSysDy(x,u,h,k,m,N,Ma,S,VtoE,maxvel,f,time,xT)
%function[rhsu]=BurgersDGrhs1D(x,u,h,k,m,N,Ma,S,VtoE,maxvel)
%Purpose:Evaluate the RHS of Burgers equations usinga DGmethod
Imat=eye(m+1);ue=zeros(N+2,2*2); %2 classes
%Extend data and assign boundary conditions
[ue]=extendDG(u(VtoE),'D',u(:,1),'D',u(:,end));
%Compute numerical fluxes at interfaces
um = ue(1:2:3,:);
up = ue(2:2:4,:);
u1 = ue(1:2,:);
u2 = ue(3:4,:);
uTm = xT(um);
uTp = xT(up);
uT = [uTm;uTp];
fluxr(1,:)=(f(u1(2,2:N+1),uTp(2:N+1),1)+f(u1(1,3:N+2),uTm(3:N+2),1))/2 ...
        -maxvel/2.*(u1(1,3:N+2)-u1(2,2:N+1));
fluxl(1,:)=(f(u1(2,1:N),uTp(1:N),1)+f(u1(1,2:N+1),uTm(2:N+1),1))/2 ...
        -maxvel/2.*(u1(1,2:N+1)-u1(2,1:N));
fluxr(2,:)=(f(u2(2,2:N+1),uTp(2:N+1),2)+f(u2(1,3:N+2),uTm(3:N+2),2))/2 ...
        -maxvel/2.*(u2(1,3:N+2)-u2(2,2:N+1));
fluxl(2,:)=(f(u2(2,1:N),uTp(1:N),2)+f(u2(1,2:N+1),uTm(2:N+1),2))/2 ...
        -maxvel/2.*(u2(1,2:N+1)-u2(2,1:N));
%Compute right hand side of Burger'sequation
ru1=S'*(f(u1(:,2:N+1),uT(:,2:N+1),1))-(Imat(:,m+1)*fluxr(1,:)-Imat(:,1)*fluxl(1,:));
ru2=S'*(f(u2(:,2:N+1),uT(:,2:N+1),2))-(Imat(:,m+1)*fluxr(2,:)-Imat(:,1)*fluxl(2,:));
rhsu(1:2,:)=(h/2*Ma)\ru1;
rhsu(3:4,:)=(h/2*Ma)\ru2;
return