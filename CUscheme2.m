function [y1,y2,dt] = CUscheme2(U1,U2,dx,q,dv,v)
    CFL = 0.5;
    n = length(U1);
    
    %Calculation
    ux1 = minmod(U1,dx);
    um1 = U1(1:end-1)+ux1(1:end-1)*dx/2; % u- j+1/2
    up1 = U1(2:end)-ux1(2:end)*dx/2; % u+ j+1/2
    
    ux2 = minmod(U2,dx);
    um2 = U2(1:end-1)+ux2(1:end-1)*dx/2; % u- j+1/2
    up2 = U2(2:end)-ux2(2:end)*dx/2; % u+ j+1/2

    um = um1+2*um2; % u- j+1/2
    up = up1+2*up2; % u+ j+1/2
    
    %vector
    bm = dv(um,1,1).*um1+v(um,1)+dv(um,2,2).*um2+v(um,2);
    cm = dv(um,1,1).*um1.*v(um,2)+dv(um,2,2).*um2.*v(um,1)+v(um,1).*v(um,2)-...
        um1.*um2.*(dv(um,2,1).*dv(um,1,2)-dv(um,2,2).*dv(um,1,1));
    eigm1 = 1/2*(bm+sqrt(bm.^2-4*cm));
    eigm2 = 1/2*(bm-sqrt(bm.^2-4*cm));
    
    bp = dv(up,1,1).*up1+v(up,1)+dv(up,2,2).*up2+v(up,2);
    cp = dv(up,1,1).*up1.*v(up,2)+dv(up,2,2).*up2.*v(up,1)+v(up,1).*v(up,2)-...
        up1.*up2.*(dv(up,2,1).*dv(up,1,2)-dv(up,2,2).*dv(up,1,1));
    eigp1 = 1/2*(bp+sqrt(bp.^2-4*cp));
    eigp2 = 1/2*(bp-sqrt(bp.^2-4*cp));
    
    a1 =[eigp1';eigm1';zeros(1,n-1)];
    a2 =[eigp2';eigm2';zeros(1,n-1)];
    ap = max(a1)'*2; 
    am = min(a2)'*2;
    
    if any(ap==am) 
        disp("ap==am")
        disp(ap(ap==am))
    end
    
    %H = (ap.*f(um)-am.*f(up)+ap.*am.*(up-um))./(ap-am);
    a = max(max(ap,abs(am)));
    H1 = (q(um1,um,1)+q(up1,up,1))/2-(a/2)*(up1-um1);%(ap.*q(um1,um,1)-am.*q(up1,up,1)+ap.*am.*(up1-um1))./(ap-am);
    H2 = (q(um2,um,2)+q(up2,up,2))/2-(a/2)*(up2-um2);%(ap.*q(um2,um,2)-am.*q(up2,up,2)+ap.*am.*(up2-um2))./(ap-am);
    
%     a = max(max(ap,abs(am)));
%     U = U1+U2*2;
%     H1 = (1/2)*(q(U1(1:end-1),U(1:end-1),1)+q(U1(2:end),U(2:end),1))...
%             -(a/2).*(U1(2:end)-U1(1:end-1));
%     H2 = (1/2)*(q(U2(1:end-1),U(1:end-1),2)+q(U2(2:end),U(2:end),2))...
%             -(a/2).*(U2(2:end)-U2(1:end-1));

    y1 = -(1/dx)*(H1(2:end)-H1(1:end-1));
    y2 = -(1/dx)*(H2(2:end)-H2(1:end-1));
    dt = CFL*dx/max(max(ap,abs(am)));
end