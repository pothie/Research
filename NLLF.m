% x: grid of distance (vector)
% t: grid of time (vector) last element:T - floor(T/dt)*dt
% u: initial condition
function U = NLLF(x,t,u,f)
    k0 = 50;
    v1max = 120; %m/s (120km/h)
    dq = @(x) v1max.*exp(-(x./k0).^2./2).*(1-x/k0);
    n_t = length(t)-1;
    %CFL = 0.9;
    dt = t(2)-t(1);
    dx = x(2)-x(1);
    U = zeros(length(x),n_t);
    U(:,1) = u;
    % alpha = max flux'
    flux = @(x,y,alpha) (1/2)*(f(x)+f(y))-(alpha/2)*(-x+y);

    for i=1:n_t
        a = max(dq(U(:,i)));
        if (i==n_t) && (t(end)-t(end-1)~=dt)
            dt = t(end)-t(end-1);
        end
        U(2:end-1,i+1)=U(2:end-1,i)-(dt/dx)*(flux(U(2:end-1,i),U(3:end,i),a)...
            -flux(U(1:end-2,i),U(2:end-1,i),a));
        
        %Periodic BC
        %U(end,i+1) = U(end,i)-(dt/dx)*(flux(U(end,i),U(2,i),dt)...
        %    -flux(U(end-1,i),U(end,i),dt));
        %U(1,i+1) = U(end,i);
        
        %One-sided BC (FTBS)
        %Ghost points
        Uend1 = 0;
        U1 = 0;
        U(end,i+1) = U(end,i)-(dt/dx)*(flux(U(end,i),Uend1,a)-flux(U(end-1,i),U(end,i),a));
        U(1,i+1) = U(1,i)-(dt/dx)*(flux(U(1,i),U(2,i),a)-flux(U1,U(1,i),a));
    end    
end