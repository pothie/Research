function[ue]=extendDG(u,BCl,ul,BCr,ur)
%Purpose:Extenddependentandindependentvectorsuwithm+1entries
%subjecttoapproproateboundaryconditionsforDGformulation
%BC="D"-Dirichlet;BC="N"-Neumann;BC="P"-periodic
%u-BCvalue-onlyactiveforDirichletBC
dim=size(u);m=dim(1)-1;N=dim(2);
ue=zeros(m+1,N+2);ue(:,2:N+1)=u;
%Periodicextensionofu
if (BCl=='P')||(BCr=='P')
ue(:,1)=u(:,N);ue(:,N+2)=u(:,1);
return;
end
%Leftextension
if BCl=='D'
ue(:,1)=ul;
else
ue(:,1)=flipud(u(:,1));
end
%Rightextension
if BCr=='D'
ue(:,N+2)=ur;%-flipud(u(:,N))+2*ur;
else
ue(:,N+2)=flipud(u(:,N));
end
return