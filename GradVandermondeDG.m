function[Vr]=GradVandermondeDG(m,r)
%function[Vr]=GradVandermondeDG(m,r)
%Purpose:InitializethegradientoftheVandermondematrix
%ofordermat(r)
Vr=zeros(length(r),(m+1));
for i=0:m
[Vr(:,i+1)]=GradLegendreP(r(:),i);
end
return