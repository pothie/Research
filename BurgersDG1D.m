function[u]=BurgersDG1D(x,u,h,m,N,CFL,FinalTime,f,df)
%function[u]=BurgersDG1D(x,u,h,m,N,CFL,FinalTime)
%Purpose:Integrate 1D equation until FinalTime using a DG
%scheme and 3rd order SSP-RKmethod
%Initialize operators at Legendre Gauss Lobatto grid
r=LegendreGL(m);V=VandermondeDG(m,r);D=DmatrixDG(m,r,V);
Ma=inv(V*V');S=Ma*D;iV=inv(V);
%Compute operator for WENO smoothness evaluator
[Q,Xm,Xp]=WENODGWeights(m,iV);
%Initialize extraction vector
VtoE=zeros(2,N);
for j=1:N
VtoE(1,j)=(j-1)*(m+1)+1;VtoE(2,j)=j*(m+1);
end
%Compute smallest spatial scale timestep
rLGLmin = min(abs(r(1)-r(2)));
time=0;tstep=0;
%integrate scheme
while(time<FinalTime)
%Decide on timestep
maxvel=max(max(abs(df(u)))); %df
k=CFL*rLGLmin*h/maxvel; 
if(time+k>FinalTime) k=FinalTime-time;end
%Update solution-stage1
rhsu=BurgersDGrhs1D(x,u,h,k,m,N,Ma,S,VtoE,maxvel,f,time);
u1=u+k*rhsu;
u1=WENOlimitDG(x,u1,m,h,N,V,iV,Q,Xm,Xp);
%Update solution-stage2
rhsu=BurgersDGrhs1D(x,u1,h,k,m,N,Ma,S,VtoE,maxvel,f,time);
u2=(3*u+u1+k*rhsu)/4;
u2=WENOlimitDG(x,u2,m,h,N,V,iV,Q,Xm,Xp);
%Update solution-stage3
rhsu=BurgersDGrhs1D(x,u2,h,k,m,N,Ma,S,VtoE,maxvel,f,time);
u=(u+2*u2+2*k*rhsu)/3;
u=WENOlimitDG(x,u,m,h,N,V,iV,Q,Xm,Xp);
% plot(x,u);
% ylim([0 80]);
% title(time);
% pause(0.1);
time=time+k;tstep=tstep+1;
end
return