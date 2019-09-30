% x: grid of distance (vector)
% t: grid of time (vector) last element:T - floor(T/dt)*dt
% u: initial condition
% f: f(u)
% df: flux
% Uend1: U(end+1,t)
% U1: U(-1,t)
function U = NLLF(x,t,u,f,df,U1,Uend1)
    if t == 0
        U = u;
    elseif x == 0
        U = U1;
    else
    n_t = length(t);
    dt = t(2)-t(1);
    dx = x(2)-x(1);
    U = zeros(length(x),n_t);
    U(:,1) = u;
    % alpha = max flux'
    flux = @(x,y,alpha) (1/2)*(f(x)+f(y))-(1*alpha)*(-x+y);

    for i=1:n_t
        a = max(abs(df(U(:,i))));
        if (i==n_t) && (t(end)-t(end-1)~=dt)
            dt = t(end)-t(end-1);
        else
        U(2:end-1,i+1)= U(2:end-1,i)-(dt/dx)*(flux(U(2:end-1,i),U(3:end,i),a)...
            -flux(U(1:end-2,i),U(2:end-1,i),a));
        
        %Periodic BC
        %U(end,i+1) = U(end,i)-(dt/dx)*(flux(U(end,i),U(2,i),dt)...
        %    -flux(U(end-1,i),U(end,i),dt));
        %U(1,i+1) = U(end,i);
        
        %One-sided BC (FTBS)
        
        %Ghost points
        if (1.125<=t(i+1)) && (t(i+1)<=1.175)
            U(end,i+1) = 200;
        else
            U(end,i+1) = ...
                U(end,i)-(dt/dx)*(flux(U(end,i),Uend1(i),a)-flux(U(end-1,i),U(end,i),a));
        end
        U(1,i+1) = U1(i); % Given U(0,t)
         %U(1,i)-(dt/dx)*(flux(U(1,i),U(2,i),a)-flux(U1(i),U(1,i),a));

        end
        
    end  
    end
end
