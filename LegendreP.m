function[P]=LegendreP(x,m);
%function[P]=LegendreP(x,m)
%Purpose:Evaluateorhonormalm'thorderLegendrePolynomialatpointx
xp=x;dims=size(xp);if(dims(2)==1)xp=xp';end;
%InitialvaluesP_0(x)andP_1(x)
PL=zeros(m+1,length(xp));
PL(1,:)=sqrt(1.0/2.0);if(m==0)P=PL';return;end;
PL(2,:)=sqrt(3.0/2.0)*xp;if(m==1)P=PL(m+1,:)';return;end;
%Forwardrecurrenceusingthesymmetryoftherecurrence.
aold=sqrt(1.0/3.0);
for i=1:m-1
    anew=2/(2*i+2)*sqrt((i+1)*(i+1)*(i+1)*(i+1)/(2*i+1)/(2*i+3));
    PL(i+2,:)=1/anew*(-aold*PL(i,:)+xp.*PL(i+1,:));
    aold=anew;
end;
P=PL(m+1,:)';
return