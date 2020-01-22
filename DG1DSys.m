function[u]=DG1DSys(x,up,ul,h,m,N,CFL,FinalTime,f,df)
%function[u]=BurgersDG1D(x,u,h,m,N,CFL,FinalTime)
%Purpose:Integrate 1D equation until FinalTime using a DG
%scheme and 3rd order SSP-RKmethod
%Initialize operators at Legendre Gauss Lobatto grid
r=LegendreGL(m);V=VandermondeDG(m,r);D=DmatrixDG(m,r,V);
Ma=inv(V*V');S=Ma*D;iV=inv(V);
%Compute operator for WENO smoothness evaluator
[Q,Xm,Xp]=WENODGWeights(m,iV);
%Initializeextractionvector
VtoE=zeros(2,N);
for j=1:N
VtoE(1,j)=(j-1)*(m+1)+1;VtoE(2,j)=j*(m+1);
end
%Initializefiltermatrix
%F=FilterDG(m,0,10,V);
%Computesmallestspatialscaletimestep
rLGLmin = min(abs(r(1)-r(2)));
time=0;tstep=0;
%Initializeparametersfornonlinearviscosity
nu=zeros(m+1,N);nu0=2;kappa=-6;c2=1;
%integratescheme
while(time<FinalTime)
%Decideontimestep
maxvel=max(max(abs(df(u))));k=CFL*rLGLmin*h/maxvel; %df
if(time+k>FinalTime) k=FinalTime-time;end
%Updatesolution-stage1
rhsu=BurgersDGrhs1D(x,u,h,k,m,N,Ma,S,VtoE,maxvel,f,time);
u1=u+k*rhsu;
u1=WENOlimitDG(x,u1,m,h,N,V,iV,Q,Xm,Xp);
%Updatesolution-stage2
rhsu=BurgersDGrhs1D(x,u1,h,k,m,N,Ma,S,VtoE,maxvel,f,time);
u2=(3*u+u1+k*rhsu)/4;
u2=WENOlimitDG(x,u2,m,h,N,V,iV,Q,Xm,Xp);
%Updatesolution-stage3
rhsu=BurgersDGrhs1D(x,u2,h,k,m,N,Ma,S,VtoE,maxvel,f,time);
u=(u+2*u2+2*k*rhsu)/3;
u=WENOlimitDG(x,u,m,h,N,V,iV,Q,Xm,Xp);
plot(x,u);
ylim([0 80]);
title(time);
pause(0.1);
time=time+k;tstep=tstep+1;
end
return