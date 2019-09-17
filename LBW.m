% x: grid of distance (vector)
% t: grid of time (vector) last element:T - floor(T/dt)*dt
% u: initial condition
function U = LBW(x,t,u)
    n_t = length(t);
    dt = t(2)-t(1);
    dx = x(2)-x(1);
    lambda = dt/dx;
    U = zeros(length(x),n_t);
    U(:,1) = u;

    for j=1:n_t
        if j == n_t
            lambda = (t(end)-t(end-1))/dx;
        end
        
        U(3:end-1,j+1)=U(3:end-1,j)-(1/2)*lambda*(3*U(3:end-1,j)-4*U(2:end-2,j)+U(1:end-3,j))...
            +(1/2)*lambda^2*(U(3:end-1,j)-2*U(2:end-2,j)+U(1:end-3,j));
        U(2,j+1)=U(2)-(1/2)*lambda*(3*U(2,j)-4*U(1,j)+U(end-1,j))...
            +(1/2)*lambda^2*(U(2,j)-2*U(1,j)+U(end-1,j));
        U(1,j+1)=U(1)-(1/2)*lambda*(3*U(1,j)-4*U(end-1,j)+U(end-2,j))...
            +(1/2)*lambda^2*(U(1,j)-2*U(end-1,j)+U(end-2,j));
        U(end,j+1)=U(1,j+1);
        
    end

end

