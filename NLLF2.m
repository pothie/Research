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
    U = zeros(length(x),ceil(T/(3e-3/60)));
    U1 = U;
    U2 = U;
    A = zeros(1,ceil(T/(3e-3/60)));
    U1(:,1) = ux0(:,1);
    U2(:,1) = ux0(:,2);
    U(:,1) = U1(:,1)+U2(:,1);
    % alpha = max eigenvalue of A
    %flux = @(x,y,alpha,n) (1/2)*(q(x,n)+q(y,n))+(alpha/2)*(x-y);
    tstep = 1;
    tgrid(tstep) = 0;
    
    while T-tgrid(tstep) > 0
        
        b = @(i) dv(U(i,tstep),1).*U1(i,tstep)+v(U(i,tstep),1)+...
            dv(U(i,tstep),2).*U2(i,tstep)+v(U(i,tstep),2);
        c = @(i) dv(U(i,tstep),1).*U1(i,tstep).*v(U(i,tstep),2)+...
            dv(U(i,tstep),2).*U2(i,tstep).*v(U(i,tstep),1)+...
            v(U(i,tstep),1).*v(U(i,tstep),2);
        eig = @(x) 1/2*(b(x)+sqrt(b(x).^2-4*c(x)));
        a = max(eig(1:length(x)));
        A(tstep) = a;
        dt = CFL*dx/a;
        %a = CFL*dt/dx;
        tpass = tgrid(tstep);
        if tpass+dt>T
            dt = T-tpass;
        end
        
        flux1 = @(x,y,step) ...
            (1/2)*(q(U1(x,step),U(x,step),1)+q(U1(y,step),U(y,step),1))+(a/2)*(U1(x,step)-U1(y,step));
        flux2 = @(x,y,step) ...
            (1/2)*(q(U2(x,step),U(x,step),2)+q(U2(y,step),U(y,step),2))+(a/2)*(U2(x,step)-U2(y,step));
        
        U1(2:end-1,tstep+1)= U1(2:end-1,tstep)-...
            (dt/dx)*(flux1(2:length(x)-1,3:length(x),tstep)...
            -flux1(1:length(x)-2,2:length(x)-1,tstep));
        U2(2:end-1,tstep+1)= U2(2:end-1,tstep)-...
            (dt/dx)*(flux2(2:length(x)-1,3:length(x),tstep)...
            -flux2(1:length(x)-2,2:length(x)-1,tstep));
        
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
        tgrid(tstep+1) = tpass+dt;
        tstep = tstep+1;
    end  
    end
end