function[rhsu]=BurgersDGrhs1D(x,u,h,k,m,N,Ma,S,VtoE,maxvel,f,time)
%function[rhsu]=BurgersDGrhs1D(x,u,h,k,m,N,Ma,S,VtoE,maxvel)
%Purpose:Evaluate the RHS of Burgers equations usinga DGmethod
Imat=eye(m+1);ue=zeros(N+2,2);
%Extend data and assign boundary conditions
[ue]=extendDG(u(VtoE),'P',0,'P',0);
%Compute numerical fluxes at interfaces
fluxr=(f(ue(2,2:N+1))+f(ue(1,3:N+2)))/2-maxvel/2.*(ue(1,3:N+2)-ue(2,2:N+1));
fluxl=(f(ue(2,1:N))+f(ue(1,2:N+1)))/2-maxvel/2.*(ue(1,2:N+1)-ue(2,1:N));
%Compute right hand side of Burger'sequation
ru=S'*(f(u))-(Imat(:,m+1)*fluxr(1,:)-Imat(:,1)*fluxl(1,:));
rhsu=(h/2*Ma)\ru;
return
