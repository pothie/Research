function dxTnew = udxT(x1,x2,xT,n,L,TH,vmax,vc,kjam,kc)
% x1,x2 are column vectors
    w = (vc*kc)/(kjam-kc);
    af = L+TH.*vmax;
    ac = TH*w*kjam;
    bf = -TH.*((vmax-vc)/kc);
    bc = L-TH*w;
    xTotal = [x1 x2];
     
    %preallocate
    dxTnew = -1*ones(size(xT));
    
    %calculation
    cindex = find(kc<=xT&xT<=kjam+1e-1); %tol 
    findex = find(0-1e-1<=xT&xT<kc); %tol
    dxTnew(findex) = (bf(n)+((af(1)-xTotal(findex,:)*bf').^2+4*bf(1)*xTotal(findex,:)*af').^(-0.5)...
        .*((xTotal(findex,:)*bf'-af(1))*bf(n)+2*bf(1)*af(n)))/(2*bf(1));
   
    dxTnew(cindex) = (bc(n)+((ac(1)-xTotal(cindex,:)*bc').^2+4*bc(1)*xTotal(cindex,:)*ac').^(-0.5)...
        .*((xTotal(cindex,:)*bc'-ac(1))*bc(n)+2*bc(1)*ac(n)))./(2*bc(1));
       
    if any(xT<0)
        disp("neg")
    end
        
end