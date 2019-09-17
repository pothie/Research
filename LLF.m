% x: grid of distance (vector)
% t: grid of time (vector) last element:T - floor(T/dt)*dt
% u: initial condition
function U = LLF(x,t,u)
    n_t = length(t);
    dt = t(2)-t(1);
    dx = x(2)-x(1);
    U = zeros(length(x),n_t);
    U(:,1) = u;

    for i=1:n_t
        if i == n_t
            dt = t(end)-t(end-1);
        end
        
        U(2:end-1,i+1)=(1/2)*(U(1:end-2,i)+U(3:end,i))...
            -(dt/(dx*2))*(-U(1:end-2,i)+U(3:end,i));
        U(end,i+1) = (1/2)*(U(2,i)+U(end-1,i))-(dt/(dx*2))*(-U(2,i)+U(end-1,i));
        U(1,i+1) = U(end,i);
    end    
end

