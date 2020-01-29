function[rhsu1,rhsu2]=DGrhs1DSys(x,u1,u2,h,k,m,N,Ma,S,VtoE,maxvel,f,time,xT)
%function[rhsu]=BurgersDGrhs1D(x,u,h,k,m,N,Ma,S,VtoE,maxvel)
%Purpose:Evaluate the RHS of Burgers equations usinga DGmethod
Imat=eye(m+1);%ue=zeros(N+2,2); 
%Extend data and assign boundary conditions
u1 = u1';
u2 = u2';
[ue1]=extendDG(u1(VtoE),'D',u1(:,1),'D',u1(:,end)); %row 
[ue2]=extendDG(u2(VtoE),'D',u2(:,1),'D',u2(:,end));
%ue1 = ue1';%column
%ue2 = ue2';
%Compute numerical fluxes at interfaces
uT = xT(ue1,ue2);
if m == 1
    fluxr(1,:)=(f(ue1(2,2:N+1),uT(2,2:N+1),1)+f(ue1(1,3:N+2),uT(1,3:N+2),1))/2 ...
            -maxvel/2.*(ue1(1,3:N+2)-ue1(2,2:N+1));
    fluxl(1,:)=(f(ue1(2,1:N),uT(2,1:N),1)+f(ue1(1,2:N+1),uT(1,2:N+1),1))/2 ...
            -maxvel/2.*(ue1(1,2:N+1)-ue1(2,1:N));
    fluxr(2,:)=(f(ue2(2,2:N+1),uT(2,2:N+1),2)+f(ue2(1,3:N+2),uT(1,3:N+2),2))/2 ...
            -maxvel/2.*(ue2(1,3:N+2)-ue2(2,2:N+1));
    fluxl(2,:)=(f(ue2(2,1:N),uT(2,1:N),2)+f(ue2(1,2:N+1),uT(1,2:N+1),2))/2 ...
            -maxvel/2.*(ue2(1,2:N+1)-ue2(2,1:N));
else
    
end
%Compute right hand side of Burger'sequation
ru1=S'*f(ue1(:,2:N+1),uT(:,2:N+1),1)-(Imat(:,m+1)*fluxr(1,:)-Imat(:,1)*fluxl(1,:));
ru2=S'*(f(ue2(:,2:N+1),uT(:,2:N+1),2))-(Imat(:,m+1)*fluxr(2,:)-Imat(:,1)*fluxl(2,:));
rhsu1=(h/2*Ma)\ru1;% row
rhsu2=(h/2*Ma)\ru2;
return