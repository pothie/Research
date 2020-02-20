% LF flux for Dynamic PCE
function [U1,U2] = LFDy(Up1,Up2,a,dx,dt,q,xT)
        Upt = xT(Up1,Up2);

        flux1 = (1/2)*(q(Up1(1:end-1),Upt(1:end-1),1)+q(Up1(2:end),Upt(2:end),1))...
            -(a/2).*(Up1(2:end)-Up1(1:end-1));
        flux2 = (1/2)*(q(Up2(1:end-1),Upt(1:end-1),2)+q(Up2(2:end),Upt(2:end),2))...
            -(a/2).*(Up2(2:end)-Up2(1:end-1));
        
        U1 =-1/dx*(flux1(2:end)-flux1(1:end-1));
        U2 =-1/dx*(flux2(2:end)-flux2(1:end-1));

end