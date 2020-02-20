function dxTnew = udxT(x1,x2,xT,n,L,TH,vmax,vc,kjam,kc)
% x1,x2 are column vectors
    w = (vc*kc)/(kjam-kc);
    [~,m] = size(x1);
    af = L+TH.*vmax; 
    af = diag(repelem(af,m));
    af = [af(1:m,1:m) af(m+1:end,m+1:end)];
    ac = TH*w*kjam; 
    ac = diag(repelem(ac,m));
    ac = [ac(1:m,1:m) ac(m+1:end,m+1:end)];
    bf = -TH.*((vmax-vc)/kc);
    bf = diag(repelem(bf,m));
    bf = [bf(1:m,1:m) bf(m+1:end,m+1:end)];
    bc = L-TH*w;
    bc = diag(repelem(bc,m));
    bc = [bc(1:m,1:m) bc(m+1:end,m+1:end)];
    xTotal = [x1 x2];
     
    %preallocate
    dxTnew = -1*ones(size(xT));
    
    %calculation
%     cindex = find(kc<=xT);%&xT<=kjam+1e-1); %tol 
%     findex = find(xT<kc); %0-1e-1<=xT&
%     dxTnew(findex) = (bf(n)+((af(1)-xTotal(findex,:)*bf').^2+4*bf(1)*xTotal(findex,:)*af').^(-0.5)...
%         .*((xTotal(findex,:)*bf'-af(1))*bf(n)+2*bf(1)*af(n)))/(2*bf(1));
%    
%     dxTnew(cindex) = (bc(n)+((ac(1)-xTotal(cindex,:)*bc').^2+4*bc(1)*xTotal(cindex,:)*ac').^(-0.5)...
%         .*((xTotal(cindex,:)*bc'-ac(1))*bc(n)+2*bc(1)*ac(n)))./(2*bc(1));
    dxTnewf = (bf(n)+((af(1)-xTotal*bf').^2+4*bf(1)*xTotal*af').^(-0.5)...
        .*((xTotal*bf'-af(1))*bf(n)+2*bf(1)*af(n)))/(2*bf(1));
   
    dxTnew = (bc(n)+((ac(1)-xTotal*bc').^2+4*bc(1)*xTotal*ac').^(-0.5)...
        .*((xTotal*bc'-ac(1))*bc(n)+2*bc(1)*ac(n)))./(2*bc(1));
       
    dxTnew(xT<kc) = dxTnewf(xT<kc); 
    
%     if any(xT<0)
%         disp("neg")
%     end
        
end