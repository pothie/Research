% x: grid of distance (vector)
% T: Total time
% ux0: initial condition [U1(0,t),U2(0,t)]
% v: speed (kTotal)
% dv: derivative of speed
% q: flow
% pce: pce values of different classes
function [U,U1,U2,tgrid] = NLLF4(x,T,ux0,v,dv,q,xT)
    CFL = 0.9;
    dx = x(2)-x(1);
    % preallocate U,U1,U2
    U = zeros(length(x),ceil(T/2));
    U1 = U;
    U2 = U;
    U1(:,1) = ux0(:,1);
    U2(:,1) = ux0(:,2);
    U(:,1) = xT(U1(:,1),U2(:,1));
    Up1 = U1(:,1);
    Up2 = U2(:,1);
    Upt = U(:,1);
    tstep = 1;
    tgrid(tstep) = 0;
    
    while T-tgrid(tstep) > 0
        b = @(i) dv(Up1(i),Up2(i),Upt(i),1,1).*Up1(i)+v(Upt(i),1)+...
            dv(Up1(i),Up2(i),Upt(i),2,2).*Up2(i)+v(Upt(i),2);
        c = @(i) dv(Up1(i),Up2(i),Upt(i),1,1).*Up1(i).*v(Upt(i),2)+...
            dv(Up1(i),Up2(i),Upt(i),2,2).*Up2(i).*v(Upt(i),1)+v(Upt(i),1).*v(Upt(i),2);%-...
            %Up1(i).*Up2(i).*(dv(Upt(i),1,1).*dv(Upt(i),2,2)-dv(Upt(i),1,2).*dv(Upt(i),2,1));
        eig = @(x) 1/2*(b(x)+sqrt(b(x).^2-4*c(x)));
        a = max(abs(eig(1:length(x))));

        flux1 = @(x,y) ...
            (1/2)*(q(Up1(x),Upt(x),1)+q(Up1(y),Upt(y),1))+(a/2)*(Up1(x)-Up1(y));
        flux2 = @(x,y) ...
            (1/2)*(q(Up2(x),Upt(x),2)+q(Up2(y),Upt(y),2))+(a/2)*(Up2(x)-Up2(y));
        
        dt = CFL*dx/a;
        tpass = tgrid(tstep);
        if tpass+dt>T
            dt = T-tpass;
        end
        
        U1(2:end-1,tstep+1)= Up1(2:end-1)-...
            (dt/dx)*(flux1(2:length(x)-1,3:length(x))...
            -flux1(1:length(x)-2,2:length(x)-1));
        U2(2:end-1,tstep+1)= Up2(2:end-1)-...
            (dt/dx)*(flux2(2:length(x)-1,3:length(x))...
            -flux2(1:length(x)-2,2:length(x)-1));
        
        U1(end,tstep+1) = U1(end-1,tstep+1);%U1(end,tstep)-...
            %(dt/dx)*(flux1(length(x),length(x))...
            %-flux1(length(x)-1,length(x)));
        U2(end,tstep+1) = U2(end-1,tstep+1);%U2(end,tstep)-...
            %(dt/dx)*(flux1(length(x),length(x))...
            %-flux1(length(x)-1,length(x)));
        
        U1(1,tstep+1) = U1(2,tstep+1);%U1(1,tstep)-...
            %(dt/dx)*(flux1(2,1) -flux1(1,1));
        U2(1,tstep+1) = U2(2,tstep+1);%U2(1,tstep)-...
            %(dt/dx)*(flux1(2,1) -flux1(1,1));
        
        U(:,tstep+1) = xT(U1(:,tstep+1),U2(:,tstep+1)); 
        
        Up1 = U1(:,tstep+1);
        Up2 = U2(:,tstep+1);
        Upt = U(:,tstep+1);
        tgrid(tstep+1) = tpass+dt;
        tstep = tstep+1;
    end  
end