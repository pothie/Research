function[u1,u2]=DG1DSys(x,u1,u2,h,m,N,CFL,FinalTime,q,v,dv,xT)
%function[u]=BurgersDG1D(x,u,h,m,N,CFL,FinalTime)
%Purpose:Integrate 1D equation until FinalTime using a DG
%scheme and 3rd order SSP-RKmethod
%Initialize operators at Legendre Gauss Lobatto grid
r=LegendreGL(m);V=VandermondeDG(m,r);D=DmatrixDG(m,r,V);
Ma=inv(V*V');S=Ma*D;
iV=inv(V);
%Compute operator for WENO smoothness evaluator
[Q,Xm,Xp]=WENODGWeights(m,iV);
%Initialize extraction vector
VtoE = reshape(1:(m+1)*N,m+1,N); 
%Compute smallest spatial scale timestep
rLGLmin = min(abs(r(1)-r(2)));
time=0;tstep=0;
%integratescheme
while(time<FinalTime)
%Decideontimestep
uT = xT(u1,u2);
b = dv(uT,1,1).*u1+v(uT,1)+dv(uT,2,2).*u2+v(uT,2);
c = dv(uT,1,1).*u1.*v(uT,2)+...
    dv(uT,2,2).*u2.*v(uT,1)+v(uT,1).*v(uT,2)-...
    u1.*u2.*(dv(uT,2,1).*dv(uT,1,2)-dv(uT,1,1).*dv(uT,2,2));
eig1 = 1/2*(b+sqrt(b.^2-4*c));
eig2 = 1/2*(b-sqrt(b.^2-4*c));
maxvel=max(max(abs([eig1 eig2])));
k=CFL*h*rLGLmin/2/maxvel;%*h/maxvel;
if(time+k>FinalTime) k=FinalTime-time;end
%Updatesolution-stage1
[rhsu1,rhsu2]=DGrhs1DSys(x,u1,u2,h,k,m,N,Ma,S,VtoE,maxvel,q,time,xT);
u1_1=u1'+k*rhsu1;%row
u2_1=u2'+k*rhsu2;
u1_1=WENOlimitDG(x,u1_1,m,h,N,V,iV,Q,Xm,Xp);
u2_1=WENOlimitDG(x,u2_1,m,h,N,V,iV,Q,Xm,Xp);
%Updatesolution-stage2
[rhsu1,rhsu2]=DGrhs1DSys(x,u1,u2,h,k,m,N,Ma,S,VtoE,maxvel,q,time,xT);
u1_2=(3*u1'+u1_1+k*rhsu1)/4;%row
u2_2=(3*u2'+u2_1+k*rhsu2)/4;
u1_2=WENOlimitDG(x,u1_2,m,h,N,V,iV,Q,Xm,Xp);
u2_2=WENOlimitDG(x,u2_2,m,h,N,V,iV,Q,Xm,Xp);
%Updatesolution-stage3
[rhsu1,rhsu2]=DGrhs1DSys(x,u1,u2,h,k,m,N,Ma,S,VtoE,maxvel,q,time,xT);
u1_3=(u1'+2*u1_2+2*k*rhsu1)/3;%row
u2_3=(u2'+2*u2_2+2*k*rhsu2)/3;
u1=WENOlimitDG(x,u1_3,m,h,N,V,iV,Q,Xm,Xp);
u2=WENOlimitDG(x,u2_3,m,h,N,V,iV,Q,Xm,Xp);
u1 = u1';
u2 = u2';
%  plot(x,u(1:2,:));
%  hold on
%  plot(x,u(3:4,:));
plot(x',xT(u1,u2)')
ylim([0 0.2]);
xlim([-6000 1000]);
title(time);
hold off
pause(0.1);
time=time+k;tstep=tstep+1;
end
return