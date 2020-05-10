function [y,dt] = CUscheme(U,dx,f,df)
    CFL = 0.5;
    n = length(U);
    %Calculation
    ux = minmod(U,dx);
    um = U(1:end-1)+ux(1:end-1)*dx/2; % u- j+1/2
    up = U(2:end)-ux(2:end)*dx/2; % u+ j+1/2
    %1D
    a = [df(up)';df(um)';zeros(1,n-1)]; 
    ap = max(a)'; % change to eigenvalues later
    am = min(a)';
    
    H = (ap.*f(um)-am.*f(up)+ap.*am.*(up-um))./(ap-am);
    H(isnan(H)) = (f(um(isnan(H)))+f(up(isnan(H))))/2;
    y = -(1/dx)*(H(2:n-1)-H(1:n-2));
    dt = CFL*dx/max(max(ap,abs(am)));
end
