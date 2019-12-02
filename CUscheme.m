function [y,dt] = CUscheme(U,dx,ux,um,up,ap,am,f,df)
    CFL = 0.5;
    n = length(U);
    %Calculation
    ux = minmod(U,dx);
    um = U(1:end-1)+ux(1:end-1)*dx/2; % u- j+1/2
    up = U(2:end)-ux(2:end)*dx/2; % u+ j+1/2
    a = [df(up)';df(um)';zeros(1,n-1)];
    ap = max(a)'; % change to eigenvalues later
    am = min(a)';
    H = @(j) (ap(j).*f(um(j))-am(j).*f(up(j))+ap(j).*am(j).*(up(j)-um(j)))./(ap(j)-am(j));

    y = -(1/dx)*(H(2:n-1)-H(1:n-2));
    dt = CFL*dx/max(max(ap,abs(am)));
end