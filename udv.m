function dvnew = udv(xT,m,vmax,vc,L1)
    %parameter    
    kjam = 1/L1;
    kc = kjam/6;
    w = (vc*kc)/(kjam-kc);

    %preacllocate
    dvnew = ones(size(xT));
    
    %assign values
    dvnew(0-1e-1<=xT&xT<=kc) =  -(vmax(m)-vc)/kc;
    dvnew(kc<xT&xT<=kjam+1e-1) =  -w.*kjam./(xT(kc<xT&xT<=kjam+1e-1).^2);%tol
   
%     if length(dvnew(0<=xT&xT<=kc))+length(dvnew(kc<xT&xT<=kjam)) ~= length(xT)
%         %disp("xT out of range in udv.")
%         disp(xT(find(dvnew==1)));
%     end
end