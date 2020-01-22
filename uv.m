function vnew = uv(xT,n,vmax,vc,L1)
    %parameter
    kjam = 1/L1;
    kc = kjam/6;
    w = (vc*kc)/(kjam-kc);
    
    %preallocate
    vnew = -1*ones(size(xT));
    
    %assign values
    xTf = xT(0-1e-1<=xT&xT<=kc);
    xTc = xT(kc<xT&xT<=kjam+1e-1); %tol
    vnew(0-1e-1<=xT&xT<=kc) = vmax(n)-((vmax(n)-vc)/kc).*xTf; 
    vnew(kc<xT&xT<=kjam+1e-1) = w.*(kjam./xTc-1); %tol
%     if length(vnew(0<=xT&xT<=kc))+length(vnew(kc<xT&xT<=kjam+1e-14)) ~= length(xT)
%         %disp("xT out of range in uv.")
%         disp(find(xT<0|xT>kjam+1e-14));
%     end
end          