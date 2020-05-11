function[ulimit]=MomentLimitDG(x,u,m,h,N,V,iV);
%functionulimit=MomentLimitDG(x,u,m,h,N,V,iV);
%Purpose:Applymomentlimitertou-anm'thorderpolynomial
eps0=1.0e-8;eps1=1.0e-8;
%Strengthofslopelimiter-Minmod:theta=1,MUSCL:theta=2
theta=2.0;
%Computecellaveragesandcellcenters
uh=iV*u;uh(2:(m+1),:)=0;uavg=V*uh;ucell=uavg(1,:);ulimit=u;
%Extendcellaverages
[ve]=extendDG(ucell,'N',0,'N',0);
%extractendvaluesandcellaveragesforeachelement
uel=u(1,:);uer=u(end,:);vj=ucell;vjm=ve(1:N);vjp=ve(3:N+2);
%Findelementsthatrequirelimiting
vel=vj-minmod1([(vj-uel)' (vj-vjm)' (vjp-vj)'])';% was xx;xx;xx
ver=vj+minmod1([(uer-vj)' (vj-vjm)' (vjp-vj)'])';% same as above
ids=(abs(vel-uel)<eps1&abs(ver-uer)<eps1);
mark=zeros(1,N);mark=(ids|mark);
%Computeexpansioncoefficients
uh=iV*u;
%Applylimitingwhenneeded
for i= m+1:-1:2
uh1=uh(i,:);uh2=uh(i-1,:);
[uh2e]=extendDG(uh2,'P',0,'P',0);uh2m=uh2e(1:N);uh2p=uh2e(3:N+2);
con=sqrt((2*i+1)*(2*i-1));
uh1=1/con*minmod1([con*uh1 theta*(uh2p-uh2) ...%was  xx;xx;xx
theta*(uh2-uh2m)]).*(1-mark)+mark.*uh1;
idsh=abs(uh1-uh(i,:))<eps0;mark=(idsh|mark);
uh(i,:)=uh1;
end
ulimit=V*uh;
return
