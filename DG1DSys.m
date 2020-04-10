function[u1,u2,tgrid]=DG1DSys(x,u1,u2,h,m,N,CFL,FinalTime,q,v,dv,xT)
%function[u]=BurgersDG1D(x,u,h,m,N,CFL,FinalTime)
%Purpose:Integrate 1D equation until FinalTime using a DG
%scheme and 3rd order SSP-RKmethod
r=LegendreGL(m);V=VandermondeDG(m,r);D=DmatrixDG(m,r,V);
Ma=inv(V*V');S=Ma*D;iV=inv(V);
%Compute operator for WENO smoothness evaluator
[Q,Xm,Xp]=WENODGWeights(m,iV);
%Initializeextractionvector
VtoE=zeros(2,N);
for j=1:N
    VtoE(1,j)=(j-1)*(m+1)+1;VtoE(2,j)=j*(m+1);
end
%Compute smallest spatial scale timestep
rLGLmin = abs(r(1)-r(2));
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
k=CFL*h*rLGLmin/maxvel;%*h/maxvel;
if(time+k>FinalTime) k=FinalTime-time;end
%Updatesolution-stage1
[rhsu1,rhsu2]=DGrhs1DSys(x,u1,u2,h,k,m,N,Ma,S,VtoE,maxvel,q,time,xT);
u1_1=u1+k*rhsu1;%row
u2_1=u2+k*rhsu2;
u1_1=WENOlimitDG(x,u1_1,m,h,N,V,iV,Q,Xm,Xp);
u2_1=WENOlimitDG(x,u2_1,m,h,N,V,iV,Q,Xm,Xp);
% u1_1 = MomentLimitDG(x,u1_1,m,h,N,V,iV);
% u2_1 = MomentLimitDG(x,u2_1,m,h,N,V,iV);
%Updatesolution-stage2
[rhsu1,rhsu2]=DGrhs1DSys(x,u1_1,u2_1,h,k,m,N,Ma,S,VtoE,maxvel,q,time,xT);
u1_2=(3*u1+u1_1+k*rhsu1)/4;%row
u2_2=(3*u2+u2_1+k*rhsu2)/4;
u1_2=WENOlimitDG(x,u1_2,m,h,N,V,iV,Q,Xm,Xp);
u2_2=WENOlimitDG(x,u2_2,m,h,N,V,iV,Q,Xm,Xp);
% u1_2 = MomentLimitDG(x,u1_2,m,h,N,V,iV);
% u2_2 = MomentLimitDG(x,u2_2,m,h,N,V,iV);
%Updatesolution-stage3
[rhsu1,rhsu2]=DGrhs1DSys(x,u1_2,u2_2,h,k,m,N,Ma,S,VtoE,maxvel,q,time,xT);
u1=(u1+2*u1_2+2*k*rhsu1)/3;%row
u2=(u2+2*u2_2+2*k*rhsu2)/3;
u1=WENOlimitDG(x,u1,m,h,N,V,iV,Q,Xm,Xp);
u2=WENOlimitDG(x,u2,m,h,N,V,iV,Q,Xm,Xp);
% u1 = MomentLimitDG(x,u1_3,m,h,N,V,iV);
% u2 = MomentLimitDG(x,u2_3,m,h,N,V,iV);


tgrid(:,tstep+1) = time;
time=time+k;tstep=tstep+1;
end
return