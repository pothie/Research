function[u]=DG1DSys(x,u,h,m,N,CFL,FinalTime,q,v,dv,xT)
%function[u]=BurgersDG1D(x,u,h,m,N,CFL,FinalTime)
%Purpose:Integrate 1D equation until FinalTime using a DG
%scheme and 3rd order SSP-RKmethod
%Initialize operators at Legendre Gauss Lobatto grid
r=LegendreGL(m);V=VandermondeDG(m,r);D=DmatrixDG(m,r,V);
Ma=inv(V*V');S=Ma*D;
iV=inv(V);
%Compute operator for WENO smoothness evaluator
[Q,Xm,Xp]=WENODGWeights(m,iV);
%Initializeextractionvector
VtoE = reshape(1:4*N,4,N); %2classes
%Computesmallestspatialscaletimestep
rLGLmin = min(abs(r(1)-r(2)));
time=0;tstep=0;
%integratescheme
while(time<FinalTime)
%Decideontimestep
um = u(1:2:3,:);
up = u(2:2:4,:);
uc1 = u(1:2,:);
uc2 = u(3:4,:);
uT(1,:) = xT(um);
uT(2,:) = xT(up);
b = dv(uT,1,1).*uc1+v(uT,1)+dv(uT,2,2).*uc2+v(uT,2);
c = dv(uT,1,1).*uc1.*v(uT,2)+...
    dv(uT,2,2).*uc2.*v(uT,1)+v(uT,1).*v(uT,2)-...
    uc1.*uc2.*(dv(uT,2,1).*dv(uT,1,2)-dv(uT,1,1).*dv(uT,2,2));
eig1 = 1/2*(b+sqrt(b.^2-4*c));
eig2 = 1/2*(b-sqrt(b.^2-4*c));
maxvel=max(max(abs([eig1 eig2])));
k=CFL*rLGLmin*h/maxvel;
if(time+k>FinalTime) k=FinalTime-time;end
%Updatesolution-stage1
rhsu=DGrhs1DSys(x,u,h,k,m,N,Ma,S,VtoE,maxvel,q,time,xT);
u1=u+k*rhsu;
u1(1:2,:)=WENOlimitDG(x,u1(1:2,:),m,h,N,V,iV,Q,Xm,Xp);
u1(3:4,:)=WENOlimitDG(x,u1(3:4,:),m,h,N,V,iV,Q,Xm,Xp);
%Updatesolution-stage2
rhsu=DGrhs1DSys(x,u1,h,k,m,N,Ma,S,VtoE,maxvel,q,time,xT);
u2=(3*u+u1+k*rhsu)/4;
u2(1:2,:)=WENOlimitDG(x,u2(1:2,:),m,h,N,V,iV,Q,Xm,Xp);
u2(3:4,:)=WENOlimitDG(x,u2(3:4,:),m,h,N,V,iV,Q,Xm,Xp);
%Updatesolution-stage3
rhsu=DGrhs1DSys(x,u2,h,k,m,N,Ma,S,VtoE,maxvel,q,time,xT);
u=(u+2*u2+2*k*rhsu)/3;
u(1:2,:)=WENOlimitDG(x,u(1:2,:),m,h,N,V,iV,Q,Xm,Xp);
u(3:4,:)=WENOlimitDG(x,u(3:4,:),m,h,N,V,iV,Q,Xm,Xp);
 plot(x,u(1:2,:));
 hold on
 plot(x,u(3:4,:));
plot(x,[xT(u(1:2:3,:));xT(u(2:2:4,:))])
ylim([0 50]);
xlim([0 2]);
title(time);
hold off
pause(0.1);
time=time+k;tstep=tstep+1;
end
return