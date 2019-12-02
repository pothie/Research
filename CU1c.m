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
    Upt = U(:,1);
    tstep = 1;
    tgrid(tstep) = 0;

    % preallocating
    ux = zeros(size(ux0));
    um = zeros(length(x)-1,1);
    up = um;
    ap = um;
    am = um;
    k1 = zeros(length(x)-2,1);
    k2 = k1;
    n = length(x);
    
    while T-tgrid(tstep) > 0
        
        %ux = minmod(Upt,dx);
        %um = Upt(1:end-1)+ux(1:end-1)*dx/2; % u- j+1/2
        %up = Upt(2:end)-ux(2:end)*dx/2; % u+ j+1/2
        %ap = max([df(up)';df(um)';ap'])'; % change to eigenvalues later
        %am = min([df(up)';df(um)';am'])';
        %H = @(j) (ap(j).*f(um(j))-am(j).*f(up(j))+ap(j).*am(j).*(up(j)-um(j)))./(ap(j)-am(j));
        
        [uav,dt] = CUscheme(Upt,dx,ux,um,up,ap,am,f,df);
        %dt = CFL*dx/max(max(ap,abs(am)));
        tpass = tgrid(tstep);
        if tpass+dt>T
            dt = T-tpass;
        end
        k1 = dt*(uav);
       
        % calculating un+k1/2
        Upt(2:end-1) = Upt(2:end-1)+k1/2;
        %ux = minmod(Upt,dx);
        %um = Upt(1:end-1)+ux(1:end-1)*dx/2; % u- j+1/2
        %up = Upt(2:end)-ux(2:end)*dx/2; % u+ j+1/2
        %dup = df(up);
        %dum = df(um);
        %ap = max([dup';dum';zeros(1,length(x)-1)])'; % change to eigenvalues later
        %am = min([dup';dum';zeros(1,length(x)-1)])';
        %H = @(j) (ap(j).*f(um(j))-am(j).*f(up(j))+ap(j).*am(j).*(up(j)-um(j)))./(ap(j)-am(j));
        [uav,~] = CUscheme(Upt,dx,ux,um,up,ap,am,f,df);
        k2 = dt*(uav);
        
        U(2:end-1,tstep+1)= U(2:end-1,tstep)+k2;
        
        U(end,tstep+1) = U(end-1,tstep+1); 
        U(1,tstep+1) = Up(tpass,pt);
        
        tgrid(tstep+1) = tpass+dt;
        tstep = tstep+1;
        Upt = U(:,tstep);
    end
    end
end