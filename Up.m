function y = Up(x)
    y = zeros(size(x));
    x0 = 0.1;
    x1 = 0.9;
    x2 = 1;
    y(x<=x0) = 400*x(x<=x0);
    y(x0<x & x<=x1) = 40;
    y(x1<x & x<=x2) = -400*x(x1<x&x<=x2)+400;
    y(x>x2) = 0;
end