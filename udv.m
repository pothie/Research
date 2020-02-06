function dvnew = udv(xT,m,vmax,vc,L1)
    %parameter    
    kjam = 1/L1;
    kc = kjam/6;
    w = (vc*kc)/(kjam-kc);

    %preacllocate
    dvnew = ones(size(xT));
    
    %assign values
%     dvnew(0-1e-1<=xT&xT<=kc) =  -(vmax(m)-vc)/kc;
%     dvnew(kc<xT&xT<=kjam+1e-1) =  -w.*kjam./(xT(kc<xT&xT<=kjam+1e-1).^2);%tol
    dvnew(xT<=kc) =  -(vmax(m)-vc)/kc;
    dvnew(kc<xT) =  -w.*kjam./(xT(kc<xT&xT<=kjam+1e-1).^2);%tol
end