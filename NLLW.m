% x: grid of distance (vector)
% t: grid of time (vector) last element:T - floor(T/dt)*dt
% u: initial condition
function U = NLLW(x,t,u,f)
    n_t = length(t)-1;
    dt = t(2)-t(1);
    dx = x(2)-x(1);
    U = zeros(length(x),n_t);
    U(:,1) = u;
    % alpha = 2u f=u^2
    flux = @(x,y,a) (1/2)*(f(x)+f(y))-(a/dx)*((x+y)/2).*(-f(x)+f(y));

    for i=1:n_t  
        if (i==n_t) && (t(end)-t(end-1)~=dt)
            dt = t(end)-t(end-1);
        end
        U(2:end-1,i+1)=U(2:end-1,i)-(dt/dx)*(flux(U(2:end-1,i),U(3:end,i),dt)...
            -flux(U(1:end-2,i),U(2:end-1,i),dt));
        
        %Periodic BC
        %U(end,i+1) = U(end,i)-(dt/dx)*(flux(U(end,i),U(2,i),dt)...
        %    -flux(U(end-1,i),U(end,i),dt));
        %U(1,i+1) = U(end,i);
        
        %One-sided BC (FTBS)
        U(end,i+1) = U(end,i)-(dt/dx)*(U(end,i)-U(end-1,i));
        U(1,i+1) = U(1,i)-(dt/dx)*(U(2,i)-U(1,i));
    end    
end

