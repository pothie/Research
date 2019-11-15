function vnew = uv(xT,n,vmax,vc,L1)
    %parameter
    kjam = 1/L1;
    kc = kjam/6;
    w = (vc*kc)/(kjam-kc);
    
    %preallocate
    vnew = -1*ones(size(xT));
    
    %assgin values
    vnew(0<=xT&xT<=kc) = vmax(n) - ((vmax(n)-vc)/kc).*xT(0<=xT&xT<=kc); 
    vnew(kc<xT&xT<=kjam+1e-14) = w.*(kjam./xT(kc<xT&xT<=kjam+1e-14)-1);
    if length(vnew(0<=xT&xT<=kc))+length(vnew(kc<xT&xT<=kjam+1e-14)) ~= length(xT)
        disp("xT out of range in uv.")
        disp(any(any(xT<0|xT>kjam+1e-14)));
    end
end          