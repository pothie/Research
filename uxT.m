function xTnew = uxT(x1,x2,L,TH,vmax,vc,kjam,kc)
 %x1,x2 are class density
    w = (vc*kc)/(kjam-kc);
    xTotal = [x1 x2];
    [m,n] = size(x1);
    af = L+TH.*vmax; 
    af = diag(repelem(af,n));
    af = [af(1:n,1:n) af(n+1:end,n+1:end)];
    ac = TH*w*kjam; 
    ac = diag(repelem(ac,n));
    ac = [ac(1:n,1:n) ac(n+1:end,n+1:end)];
    bf = -TH.*((vmax-vc)/kc);
    bf = diag(repelem(bf,n));
    bf = [bf(1:n,1:n) bf(n+1:end,n+1:end)];
    bc = L-TH*w;
    bc = diag(repelem(bc,n));
    bc = [bc(1:n,1:n) bc(n+1:end,n+1:end)];
    
    pf = -1*ones(m,n);
    pc = pf;
    pf = (af(1)-xTotal*bf'-((af(1)-xTotal*bf').^2+4*bf(1)*xTotal*af').^0.5)./(-2*bf(1));
    pc = (ac(1)-xTotal*bc'-((ac(1)-xTotal*bc').^2+4*bc(1)*xTotal*ac').^0.5)./(-2*bc(1));
    pf(abs(imag(pf))>1e-15) = -2; %drop complex solutions 
    
    % solution always exists during congestion
    xTnew = pc;
    %xTnew(kc<=pc) = pc(kc<=pc);
    xTnew(-1e-1<pf&pf<kc) = pf(-1e-1<pf&pf<kc);
end