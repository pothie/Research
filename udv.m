function dvnew = udv(xT,m,vmax,vc,L1)
    %parameter    
    kjam = 1/L1;
    kc = kjam/6;
    w = (vc*kc)/(kjam-kc);

    %preacllocate
    dvnew = ones(size(xT));
    
    %assign values
    dvnew(0<=xT&xT<=kc) =  -(vmax(m)-vc)/kc;
    dvnew(kc<xT&xT<=kjam+1e-14) =  -w.*(kjam./(xT(kc<xT&xT<=kjam+1e-14).^2));
    %dvnew = dvnew.*dxT(x1,x2,xT,n);
    if length(dvnew(0<=xT&xT<=kc))+length(dvnew(kc<xT&xT<=kjam+1e-14)) ~= length(xT)
        disp("xT out of range in udv.")
        disp(find(dvnew==1))
    end
end