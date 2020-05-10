% x: grid of distance (vector)
% t: grid of time (vector) last element:T - floor(T/dt)*dt
% u: initial condition
% f: f(u)
% df: f'(u)
function [U, tgrid]= CU(x,T,ux0,f,df)
    pt = [0.25 10;0.5 50;1 50;1.25 10];
    if T == 0
        U = ux0;
    elseif x == 0
        U = Up(0:1.5e-3:T,pt);
    else
    s = @(x,t) -sin(x).*sin(t)+sin(x).*cos(x).*cos(t).^2;
    CFL = 0.1;
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
    k1 = um;
    k2 = um;
    while T-tgrid(tstep) > 0
        dt = CFL*dx;
        tpass = tgrid(tstep);
        if tpass+dt>T
            dt = T-tpass;
        end
        ux = minmod(Upt,dx);
        um = Upt(1:end-1)+ux(1:end-1)*dx/2; % u- j+1/2
        up = Upt(2:end)-ux(2:end)*dx/2; % u+ j+1/2
        ap = max([df(up)';df(um)';ap'])'; % change to eigenvalues later
        am = min([df(up)';df(um)';am'])';
        H = @(j) (ap(j).*f(um(j))-am(j).*f(up(j))+ap(j).*am(j).*(up(j)-um(j)))./(ap(j)-am(j));
        k1(1) = dt*(s(x(1),tpass)-(1/(dx))*(H(1)-H(length(x)-1)));
        k1(2:end) = dt*(s(x(2:end-1),tpass)-(1/(dx))*(H(2:length(x)-1)-H(1:length(x)-2)));
        
        % calculating un+k1/2
        Upt(1:end-1) = Upt(1:end-1)+k1/2;
        ux = minmod(Upt,dx);
        um = Upt(1:end-1)+ux(1:end-1)*dx/2; % u- j+1/2
        up = Upt(2:end)-ux(2:end)*dx/2; % u+ j+1/2
        ap = max([df(up)';df(um)';zeros(1,length(x)-1)])'; % change to eigenvalues later
        am = min([df(up)';df(um)';zeros(1,length(x)-1)])';
        H = @(j) (ap(j).*f(um(j))-am(j).*f(up(j))+ap(j).*am(j).*(up(j)-um(j)))./(ap(j)-am(j));
        k2(1) = dt*(s(x(1),tpass+dt/2)-(1/(dx))*(H(1)-H(length(x)-1)));
        k2(2:end) = dt*(s(x(2:end-1),tpass+dt/2)-(1/(dx))*(H(2:length(x)-1)-H(1:length(x)-2)));
        
        U(1:end-1,tstep+1)= U(1:end-1,tstep)+k2;
        
        U(end,tstep+1) = U(1,tstep+1); 
         
        tgrid(tstep+1) = tpass+dt;
        tstep = tstep+1;
        Upt = U(:,tstep);
    end
    end
end
