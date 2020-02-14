function[rhsu1,rhsu2]=DGrhs1DSys(x,u1,u2,h,k,m,N,Ma,S,VtoE,maxvel,f,time,xT)
%function[rhsu]=BurgersDGrhs1D(x,u,h,k,m,N,Ma,S,VtoE,maxvel)
%Purpose:Evaluate the RHS of Burgers equations usinga DGmethod
Imat=eye(m+1);%ue=zeros(N+2,2); 
%Extend data and assign boundary conditions
u1 = u1'; %row
u2 = u2';
%[ue1]=extendDG(u1(VtoE),'D',u1(:,1),'D',u1(:,end)); %row 
%[ue2]=extendDG(u2(VtoE),'D',u2(:,1),'D',u2(:,end));
%ue1 = ue1';%column
%ue2 = ue2';
%Compute numerical fluxes at interfaces
% uT = xT(ue1,ue2);
uT = xT(u1,u2);
% fluxr1 = zeros(1,length(u1));
% fluxl1 = zeros(size(fluxr1));
% fluxr2 = zeros(size(fluxr1));
% fluxl1 = zeros(size(fluxr));
    fluxr = zeros(2,length(u1));
    fluxl = zeros(size(fluxr));
    fluxr(1,2:end-1)=(f(u1(end,2:end-1),uT(end,2:end-1),1)+f(u1(1,3:end),uT(1,3:end),1))/2 ...
            -maxvel/2.*(u1(1,3:end)-u1(end,2:end-1));
    fluxl(1,2:end-1)=(f(u1(end,1:end-2),uT(end,1:end-2),1)+f(u1(1,2:end-1),uT(1,2:end-1),1))/2 ...
            -maxvel/2.*(u1(1,2:end-1)-u1(end,1:end-2));
    %flux1 = (f(u1(end,1:end-1),uT(end,1:end-1),1)+f(u1(1,2:end),uT(1,2:end),1))/2 ...
            -maxvel/2.*(u1(1,2:end-1)-u1(end,1:end-2));
        
        
    fluxr(2,2:end-1)=(f(u2(end,2:end-1),uT(end,2:end-1),2)+f(u2(1,3:end),uT(1,3:end),2))/2 ...
            -maxvel/2.*(u2(1,3:end)-u2(end,2:end-1));
    fluxl(2,2:end-1)=(f(u2(end,1:end-2),uT(end,1:end-2),2)+f(u2(1,2:end-1),uT(1,2:end-1),2))/2 ...
            -maxvel/2.*(u2(1,2:end-1)-u2(end,1:end-2));
        
    fluxr(:,1) = fluxl(:,2);
    fluxl(:,end) = fluxr(:,end-1);
    
    fluxl(:,1) = fluxl(:,2);%f(u1(2,1),uT(2,1),1);
    %fluxl(:,1) = f(u1(2,1),uT(2,1),1);
    fluxr(:,end) = fluxr(:,end-1);%f(u2(2,end),uT(2,end),2);

ru1=S'*f(u1,uT,1)-(Imat(:,m+1)*fluxr(1,:)-Imat(:,1)*fluxl(1,:));
ru2=S'*f(u2,uT,2)-(Imat(:,m+1)*fluxr(2,:)-Imat(:,1)*fluxl(2,:));
rhsu1=(h/2*Ma)\ru1;% row
rhsu2=(h/2*Ma)\ru2;
return