function vnew = uv(xT,n,vmax,vc,L1)
    %parameter
    kjam = 1/L1;
    kc = kjam/6;
    w = (vc*kc)/(kjam-kc);
    
    %preallocate
    vnew = -1*ones(size(xT));
    
    %assign values    
    xTf = xT(xT<=kc);
    xTc = xT(kc<xT); 
    vnew(xT<=kc) = vmax(n)-((vmax(n)-vc)/kc).*xTf; %0-1e-1<=xT&
    vnew(kc<xT) = w.*(kjam./xTc-1); %&xT<=kjam+1e-1
end          