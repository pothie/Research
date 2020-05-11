% x: grid of distance (vector)
% t: grid of time (vector) last element:T - floor(T/dt)*dt
% u: initial condition
% f: f(u)
% df: f'(u)
function [U, tgrid]= NLLF(x,T,u,f,df)
    pt = [0 0;0.1 40;0.9 40;1 0];
    if T == 0
        U = u;
    elseif x == 0
        U = Up(0:1.5e-3:T,pt);
    else
    %s = @(x,t) -sin(x).*sin(t)+sin(x).*cos(x).*cos(t).^2;
    CFL = 0.9;
    dx = x(2)-x(1);
    U(:,1) = u;
    flux = @(x,y,alpha) (1/2)*(f(x)+f(y))+(alpha/2)*(x-y);
    tstep = 1;
    tgrid(tstep) = 0;
    
    while T-tgrid(tstep) > 0
        a = max(abs(df(U(:,tstep))));
        dt = CFL*dx/a;
        tpass = tgrid(tstep);
        if tpass+dt>T
            dt = T-tpass;
        end
        
        U(2:end-1,tstep+1)= U(2:end-1,tstep)-...
            (dt/dx)*(flux(U(2:end-1,tstep),U(3:end,tstep),a)-...
            flux(U(1:end-2,tstep),U(2:end-1,tstep),a));
        
        %U1 = Up(tpass,pt);
%         if (1.125<=tpass+dt) && (tpass+dt<=1.175)
%            Uend1 = 200;
%         elseif (tpass+dt>=1.175)
%            Uend1 = 0;
%         else 
%            Uend1 = U(end-1,tstep+1);
%         end
        
        U(end,tstep+1) = U(end-1,tstep+1);
        
        U(1,tstep+1) = Up(tpass,pt);

        tgrid(tstep+1) = tpass+dt;
        tstep = tstep+1;
    end  
    end
end
