function[V]=VandermondeDG(m,r)
%function[V]=VandermondeDG(m,r)
%Purpose:Initializethe1DVandermondeMatrix,V_{ij}=phi_j(r_i);
V=zeros(length(r),m+1);
for j=1:m+1
V(:,j)=LegendreP(r(:),j-1);
end;
return