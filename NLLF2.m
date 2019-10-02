% x: grid of distance (vector)
% t: grid of time (vector) last element:T - floor(T/dt)*dt
% u: initial condition
% f: f(u)
% df: flux
% Uend1: U(end+1,t)
% U1: U(-1,t)
function [U,U1,U2,tgrid] = NLLF2(x,T,ux0,v,dv,q,ut0)
    pt = [0.25 10;0.5 50;1 50;1.25 10];
    if T == 0
        U = ut0;
    elseif x == 0
        U1(:,1) = 0.2*Up(0:1.5e-3:T,pt);
        U2(:,1) = 0.8*Up(0:1.5e-3:T,pt);
        U = U1+U2;
    else
    CFL = 0.9;
    dx = x(2)-x(1);
    U1(:,1) = ux0(:,1);
    U2(:,1) = ux0(:,2);
    U(:,1) = U1(:,1)+U2(:,1);
    % alpha = max eigenvalue of A
    flux = @(x,y,alpha,n) (1/2)*(q(x,n)+q(y,n))+(alpha/2)*(x-y);
    tstep = 1;
    tgrid(tstep) = 0;
    
    while T-tgrid(tstep) > 0
        %dv1dk = max(dv(U(:,tstep),1));
        %dv2dk = max(dv(U(:,tstep),2));
        %eigv = 0;
        for  i=1:length(x)
            A = [dv(U(i,tstep),1).*U1(i,tstep)+v(U(i,tstep),1) dv(U(i,tstep),1).*U1(i,tstep);
                 dv(U(i,tstep),2).*U2(i,tstep) dv(U(i,tstep),2).*U2(i,tstep)+v(U(i,tstep),2)];
            [~,D] = eig(A);
            eigv = max(abs(diag(D)));
        end
        a = max(abs(eigv));
        dt = CFL*dx/a;
        %a = CFL*dt/dx;
        tpass = tgrid(tstep);
        if tpass+dt>T
            dt = T-tpass;
        end
        
        U1(2:end-1,tstep+1)= U1(2:end-1,tstep)-...
            (dt/dx)*(flux(U1(2:end-1,tstep),U1(3:end,tstep),a,1)...
            -flux(U1(1:end-2,tstep),U1(2:end-1,tstep),a,1));
        U2(2:end-1,tstep+1)= U2(2:end-1,tstep)-...
            (dt/dx)*(flux(U2(2:end-1,tstep),U2(3:end,tstep),a,2)...
            -flux(U2(1:end-2,tstep),U2(2:end-1,tstep),a,2));
        
        %If an accident happens between t=1.125 and t=1.175
        Uend1 = [0.8 0.2].*Up(tpass,pt);
        Un1 = Uend1;
        %if (1.125<=tpass+dt) && (tpass+dt<=1.175)
        %     U1(end,tstep+1) = 0.8*200;
        %     U2(end,tstep+1) = 0.2*200;
        %else
            U1(end,tstep+1) = ...
                U1(end,tstep)-(dt/dx)*(flux(U1(end,tstep),Uend1(1),a,1)...
                -flux(U1(end-1,tstep),U1(end,tstep),a,1));
            U2(end,tstep+1) = ...
                U2(end,tstep)-(dt/dx)*(flux(U2(end,tstep),Uend1(2),a,2)...
                -flux(U2(end-1,tstep),U2(end,tstep),a,2));
        %end
        U1(1,tstep+1) = Un1(1); % Given U(0,t)
        U2(1,tstep+1) = Un1(2);
                %U(1,tstep)-(dt/dx)*(flux(U(1,tstep),U(2,tstep),a)...
                %-flux(U1,U(1,tstep),a));
        U(:,tstep+1) = U1(:,tstep+1)+U2(:,tstep+1); 
        tgrid(tstep+1) = tpass+dt;
        tstep = tstep+1;
    end  
    end
end