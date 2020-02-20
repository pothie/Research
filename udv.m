function dvnew = udv(xT,m,vmax,vc,L1)
    %parameter    
    kjam = 1/L1;
    kc = kjam/6;
    w = (vc*kc)/(kjam-kc);

    %preacllocate
    dvnew = ones(size(xT));
    
    %assign values
    dvnew(xT<=kc) =  -(vmax(m)-vc)/kc;
    dvnew(kc<xT) =  -w.*kjam./(xT(kc<xT).^2);
end