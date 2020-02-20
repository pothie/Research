% x: grid of distance (vector)
% T: Total time
% ux0: initial condition [U1(0,t),U2(0,t)]
% v: speed (kTotal)
% dv: derivative of speed
% q: flow
% pce: pce values of different classes
function [U,U1,U2,tgrid] = NLLF3(x,T,ux0,v,dv,q,pce)
    %pt = [0.25 10;0.5 50;1 50;1.25 10];
    CFL = 0.9;%0.5
    dx = x(2)-x(1);
    % preallocate U,U1,U2
    U = zeros(length(x),1);
    U1 = U;
    U2 = U;
    U1(:,1) = ux0(:,1);
    U2(:,1) = ux0(:,2);
    U(:,1) = pce*[U1(:,1) U2(:,1)]';
    Up1 = U1(:,1);
    Up2 = U2(:,1);
    Upt = U(:,1);
    tstep = 1;
    tgrid(tstep) = 0;
    
    while T-tgrid(tstep) > 0
        b = dv(Upt,1,1).*Up1+v(Upt,1)+dv(Upt,2,2).*Up2+v(Upt,2);
        c = dv(Upt,1,1).*Up1.*v(Upt,2)+...
            dv(Upt,2,2).*Up2.*v(Upt,1)+v(Upt,1).*v(Upt,2)-...
            Up1.*Up2.*(dv(Upt,2,1).*dv(Upt,1,2)-dv(Upt,1,1).*dv(Upt,2,2));
        eig(1,:) = 1/2*(b+sqrt(b.^2-4*c));
        eig(2,:) = 1/2*(b-sqrt(b.^2-4*c));
        a =  max(max(abs(eig)));

        dt = CFL*dx/a;%
        tpass = tgrid(tstep);
        if tpass+dt>T
            dt = T-tpass;
        end
        
        flux1 = (1/2)*(q(Up1(1:end-1),Upt(1:end-1),1)+q(Up1(2:end),Upt(2:end),1))...
            -(a/2).*(Up1(2:end)-Up1(1:end-1));
        flux2 = (1/2)*(q(Up2(1:end-1),Upt(1:end-1),2)+q(Up2(2:end),Upt(2:end),2))...
            -(a/2).*(Up2(2:end)-Up2(1:end-1));

        U1(2:end-1,tstep+1)= Up1(2:end-1)-...
            (dt/dx)*(flux1(2:end)-flux1(1:end-1));
        U2(2:end-1,tstep+1)= Up2(2:end-1)-...
            (dt/dx)*(flux2(2:end)-flux2(1:end-1));
        
        U1(end,tstep+1) = U1(end,tstep)-...
            (dt/dx)*(q(Up1(end),Upt(end),1)-q(Up1(end-1),Upt(end-1),1));
        U2(end,tstep+1) = U2(end,tstep)-...
            (dt/dx)*(q(Up2(end),Upt(end),2)-q(Up2(end-1),Upt(end-1),2));
        
        U1(1,tstep+1) = U1(1,tstep)-...
            (dt/dx)*(q(Up1(2),Upt(2),1)-q(Up1(1),Upt(1),1));%quantity supposed to be 0
        U2(1,tstep+1) = U2(1,tstep)-...
            (dt/dx)*(q(Up2(2),Upt(2),2)-q(Up2(1),Upt(1),2));
        
        U(:,tstep+1) = pce*[U1(:,tstep+1) U2(:,tstep+1)]'; 
        
        Up1 = U1(:,tstep+1);
        Up2 = U2(:,tstep+1);
        Upt = U(:,tstep+1);
        tgrid(tstep+1) = tpass+dt;
        tstep = tstep+1;
    end  
end