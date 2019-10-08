% x: grid of distance (vector)
% t: grid of time (vector) last element:T - floor(T/dt)*dt
% u: initial condition
% f: f(u)
% df: f'(u)
function [U, tgrid]= NLLF(x,T,u,f,df)
    pt = [0.25 10;0.5 50;1 50;1.25 10];
    if T == 0
        U = u;
    elseif x == 0
        U = Up(0:1.5e-3:T,pt);
    else
    CFL = 0.9;
    dx = x(2)-x(1);
    U = zeros(length(x),ceil(T/(1.5e-3/60)));
    U(:,1) = u;
    % alpha = max q'
    flux = @(x,y,alpha) (1/2)*(f(x)+f(y))+(alpha/2)*(x-y);
    tstep = 1;
    tgrid(tstep) = 0;
    
    while T-tgrid(tstep) > 0
        a = max(abs(df(U(:,tstep))));
        dt = CFL*dx/a;
        %a = 0.9*dx/dt;
        tpass = tgrid(tstep);
        if tpass+dt>T
            dt = T-tpass;
        end
        
        U(2:end-1,tstep+1)= U(2:end-1,tstep)-...
            (dt/dx)*(flux(U(2:end-1,tstep),U(3:end,tstep),a)...
            -flux(U(1:end-2,tstep),U(2:end-1,tstep),a));
        
        U1 = Up(tpass,pt);
        %If an accident happens between t=1.125 and t=1.175
        if (1.125<=tpass+dt) && (tpass+dt<=1.175)
            Uend1 = 200;
        elseif (tpass+dt>=1.175)
            Uend1 = 0;
        else 
            Uend1 = U(end,tstep);
        end
        
        U(end,tstep+1) = Uend1;%...
                %U(end,tstep)-(dt/dx)*(flux(U(end,tstep),Uend1,a)...
                %-flux(U(end-1,tstep),U(end,tstep),a));
        
        U(1,tstep+1) = U1; % Given U(0,t)
                %U(1,tstep)-(dt/dx)*(flux(U(1,tstep),U(2,tstep),a)...
                %-flux(U1,U(1,tstep),a));
         
        tgrid(tstep+1) = tpass+dt;
        tstep = tstep+1
        
    end  
    end
end