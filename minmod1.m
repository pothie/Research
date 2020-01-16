function psi=minmod1(v)
%functionpsi=minmod(v)
%Purpose:Implementthemidmodfunctionvisavector
N=size(v,1);m=size(v,2);psi=zeros(N,1);
s=sum(sign(v),2)/m;ids=find(abs(s)==1);
if(~isempty(ids))
psi(ids)=s(ids).*min(abs(v(ids,:)),[],2);
end
return;