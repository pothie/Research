% u: column vector
% dx: scalar: delta x
function ux = minmod(u,dx)
    ux = zeros(size(u));
    a = zeros(3,1);
    for j = 2:length(u)-1
        a(1) = 2*(u(j)-u(j-1))/dx;
        a(2) = (u(j+1)-u(j-1))/(2*dx);
        a(3) = 2*(u(j+1)-u(j))/dx;
        if all(a>0)==1
            ux(j) = min(a);
        elseif all(a<0)==1
            ux(j) = max(a);
        end
    end
    ux(1) = (u(2)-u(1))/dx;
    ux(end) = (u(end)-u(end-1))/dx;
end