% x: grid of distance (vector)
% t: grid of time (vector) last element:T - floor(T/dt)*dt
% u: initial condition
% f: f(u)
% df: f'(u)
function [U, tgrid]= CU1c(x,T,ux0,f,df)
    pt = [0.25 10;0.5 50;1 50;1.25 10];
    if T == 0
        U = ux0;
    elseif x == 0
        U = Up(0:1.5e-3:T,pt);
    else
    %CFL = 0.5;
    dx = x(2)-x(1);
    U(:,1) = ux0;
    tstep = 1;
    tgrid(tstep) = 0;
    k1 = zeros(length(x)-2,1);
    k2 = k1;
    
    while T-tgrid(tstep) > 0
        Upt = U(:,tstep);
        Uptc = Upt;
        [uav,dt] = CUscheme(Upt,dx,f,df);
        tpass = tgrid(tstep);
        if tpass+dt>T
            dt = T-tpass;
        end
        % RK2 k1
        k1 = dt*(uav);
       
        % calculating k2 = un+k1/2
        Upt(2:end-1) = Upt(2:end-1)+k1/2;
        [uav,~] = CUscheme(Upt,dx,f,df);
        k2 = dt*(uav);
        U(2:end-1,tstep+1)= Uptc(2:end-1)+k2;
        
        % Boundary points U(1,:) = Up(tpass,pt)
        U(end,tstep+1) = U(end-1,tstep+1); 
        U(1,tstep+1) = Up(tpass,pt);%U(2,tstep+1);%Up(tpass,pt);
        
        tgrid(tstep+1) = tpass+dt;
        tstep = tstep+1;
    end
    end
end