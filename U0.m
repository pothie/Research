function y = U0(x,x0,ul,ur)
    y = zeros(size(x));
    y(x<=x0) = ul;
    y(x>x0) = ur;
end