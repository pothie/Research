% x: grid of distance (vector)
% T: Total time
% ux0: initial condition [U1(0,t),U2(0,t)]
% v: speed (kTotal)
% dv: derivative of speed
% q: flow
% pce: pce values of different classes
function [U,U1,U2,tgrid] = NLLF4(x,T,ux0,v,dv,q,xT)
    CFL = 0.1;
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
        b = dv(Up1,Up2,Upt,1,1).*Up1+v(Upt,1)+dv(Up1,Up2,Upt,2,2).*Up2+v(Upt,2);
        c = dv(Up1,Up2,Upt,1,1).*Up1.*v(Upt,2)+...
            dv(Up1,Up2,Upt,2,2).*Up2.*v(Upt,1)+v(Upt,1).*v(Upt,2)-...
            Up1.*Up2.*(dv(Up1,Up2,Upt,2,1).*dv(Up1,Up2,Upt,1,2)-dv(Up1,Up2,Upt,1,1).*dv(Up1,Up2,Upt,2,2));
        eig = 1/2*(b+sqrt(b.^2-4*c));
        a = max(abs(eig));
        dt = CFL*dx/a;
        
        tpass = tgrid(tstep);
        if tpass+dt>T
            dt = T-tpass;
        end
        
        [uav1,uav2] = LFDy(Up1,Up2,a,dx,dt,q,xT);
        
        tpass = tgrid(tstep);
        if tpass+dt>T
            dt = T-tpass;
        end
        % RK2 k1
        k1_1 = dt*(uav1);
        k1_2 = dt*(uav2);
       
        % calculating k2 = un+k1/2
        Up1(2:end-1) = Up1(2:end-1)+k1_1;
        Up2(2:end-1) = Up2(2:end-1)+k1_2;
        Up1(end) = Up1(end-1);
        Up2(end) = Up2(end-1);
        Up1(1) = Up1(2);
        Up2(1) = Up2(2);
        
        [uav1,uav2] = LFDy(Up1,Up2,a,dx,dt,q,xT);
        k2_1 = dt*(uav1);
        k2_2 = dt*(uav2);
        Up1(2:end-1)= U1(2:end-1,tstep)*0.75+Up1(2:end-1)*0.25+k2_1/4;
        Up2(2:end-1)= U2(2:end-1,tstep)*0.75+Up2(2:end-1)*0.25+k2_2/4;
        Up1(end) = Up1(end-1);
        Up2(end) = Up2(end-1);
        Up1(1) = Up1(2);
        Up2(1) = Up2(2);
        
        %put BC here, missing?
        [uav1,uav2] = LFDy(Up1,Up2,a,dx,dt,q,xT);
        k3_1 = dt*(uav1);
        k3_2 = dt*(uav2);

        Up1(2:end-1)= U1(2:end-1,tstep)/3+Up1(2:end-1)*2/3+k3_1*2/3;
        Up2(2:end-1)= U2(2:end-1,tstep)/3+Up2(2:end-1)*2/3+k3_2*2/3;

        U1(2:end-1,tstep+1) = Up1(2:end-1);
        U2(2:end-1,tstep+1) = Up2(2:end-1);
        Up1(end) = Up1(end-1);
        Up2(end) = Up2(end-1);
        Up1(1) = Up1(2);
        Up2(1) = Up2(2);
        
        %Boundary points 
        U1(end,tstep+1) = Up1(end-1);%Up1(end);
        U2(end,tstep+1) = Up2(end-1);%Up2(end);

        U1(1,tstep+1) = Up1(2);%Up1(1);
        U2(1,tstep+1) = Up2(2);%Up2(1);
        
        U(:,tstep+1) = xT(U1(:,tstep+1),U2(:,tstep+1)); 
        
        if mod(tpass,10)<dt
            plot(x,U(:,tstep+1));
            title(tpass);
            pause(0.1);
        end
        
        Up1 = U1(:,tstep+1);
        Up2 = U2(:,tstep+1);
        Upt = U(:,tstep+1);
        tgrid(tstep+1) = tpass+dt;
        tstep = tstep+1;
    end  
end