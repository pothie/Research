% x: x-values
% pt: 4x2 matrix with 4 turning points
function y = Up(x,pt)
    y = zeros(size(x));
    x1 = pt(1,1);
    x2 = pt(2,1);
    x3 = pt(3,1);
    x4 = pt(4,1);
    y1 = pt(1,2);
    y2 = pt(2,2);
    y(x<=x1) = y1;
    y(x1<x & x<=x2) = (y2-y1)/(x2-x1)*(x(x1<x & x<=x2)-x1)+y1;
    y(x2<x & x<=x3) = y2;
    y(x3<x & x<=x4) = (y1-y2)/(x4-x3)*(x(x3<x & x<=x4)-x3)+y2;
    y(x>x4) = y1;
end
