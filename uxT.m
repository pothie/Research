function xTnew = uxT(x1,x2,L,TH,vmax,vc,kjam,kc)
 %x1,x2 are strictly column vectors
    w = (vc*kc)/(kjam-kc);
    af = L+TH.*vmax;
    ac = TH*w*kjam;
    bf = -TH.*((vmax-vc)/kc);
    bc = L-TH*w;
    xTotal = [x1 x2];
    
    pf = -1*ones(size(x1));
    pc = -1*ones(size(x1));
    pf = (af(1)-xTotal*bf'-((af(1)-xTotal*bf').^2+4*bf(1)*xTotal*af').^0.5)./(-2*bf(1));
    pc = (ac(1)-xTotal*bc'-((ac(1)-xTotal*bc').^2+4*bc(1)*xTotal*ac').^0.5)./(-2*bc(1));
    pf(abs(imag(pf))>1e-16) = -2; %drop complex solutions
    
    % solution always exists during congestion
    xTnew = pc;%-1*ones(size(x1));
    xTnew(kc<=pc&pc<=kjam+1e-14) = pc(kc<=pc&pc<=kjam+1e-14);
    xTnew(0-1e-1<=pf&pf<kc) = pf(0-1e-1<=pf&pf<kc);
    
    % Correction
     %xTnew(xTnew<0)= 0;
     %xTnew(xTnew>kjam)= kjam;
%     if  any(xTnew>kjam+1e-14) || any(xTnew<0-1e-14)
%         disp("out of range")
%     end
end