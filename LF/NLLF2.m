% x: grid of distance (vector)
% t: grid of time (vector) last element:T - floor(T/dt)*dt
% u: initial condition
% f: f(u)
% df: flux
% Uend1: U(end+1,t)
% U1: U(-1,t)
function [U,U1,U2,tgrid] = NLLF2(x,T,ux0,v,dv,q)
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
    U = zeros(length(x),ceil(T/(3e-3/60)));
    U1 = U;
    U2 = U;
    U1(:,1) = ux0(:,1);
    U2(:,1) = ux0(:,2);
    U(:,1) = U1(:,1)+U2(:,1);
    Up1 = U1(:,1);
    Up2 = U2(:,1);
    Upt = U(:,1);
    tstep = 1;
    tgrid(tstep) = 0;
    
    while T-tgrid(tstep) > 0
        
        b = @(i) dv(Upt(i),1).*Up1(i)+v(Upt(i),1)+...
            dv(Upt(i),2).*Up2(i)+v(Upt(i),2);
        c = @(i) dv(Upt(i),1).*Up1(i).*v(Upt(i),2)+...
            dv(Upt(i),2).*Up2(i).*v(Upt(i),1)+v(Upt(i),1).*v(Upt(i),2);
        eig = @(x) 1/2*(b(x)+sqrt(b(x).^2-4*c(x)));
        a = max(eig(1:length(x)));
        dt = CFL*dx/a;
        tpass = tgrid(tstep);
        if tpass+dt>T
            dt = T-tpass;
        end
        
        flux1 = @(x,y) ...
            (1/2)*(q(Up1(x),Upt(x),1)+q(Up1(y),Upt(y),1))+(a/2)*(Up1(x)-Up1(y));
        flux2 = @(x,y) ...
            (1/2)*(q(Up2(x),Upt(x),2)+q(Up2(y),Upt(y),2))+(a/2)*(Up2(x)-Up2(y));
        
        U1(2:end-1,tstep+1)= Up1(2:end-1)-...
            (dt/dx)*(flux1(2:length(x)-1,3:length(x))...
            -flux1(1:length(x)-2,2:length(x)-1));
        U2(2:end-1,tstep+1)= Up2(2:end-1)-...
            (dt/dx)*(flux2(2:length(x)-1,3:length(x))...
            -flux2(1:length(x)-2,2:length(x)-1));
        
        %If an accident happens between t=1.125 and t=1.175
        
        if (1.125<=tpass+dt) && (tpass+dt<=1.175)
            Uend1 = 200;
        elseif (tpass+dt>=1.175)
            Uend1 = 0;
        else 
            Uend1 = U(end,tstep);
        end
        
        U1(end,tstep+1) = 0.8*Uend1;
        U2(end,tstep+1) = 0.2*Uend1;
        
        Un1 = [0.8 0.2].*Up(tpass,pt);
        U1(1,tstep+1) = Un1(1); % Given U(0,t)
        U2(1,tstep+1) = Un1(2);
        
        U(:,tstep+1) = U1(:,tstep+1)+U2(:,tstep+1); 
        
        Up1 = U1(:,tstep+1);
        Up2 = U2(:,tstep+1);
        Upt = U(:,tstep+1);
        tgrid(tstep+1) = tpass+dt;
        tstep = tstep+1;
    end  
    end
end
