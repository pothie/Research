% x: x-values
% pt: mx2 matrix with m turning points
function y = Up(x,pt)
    y = zeros(size(x));
    [m,~] = size(pt);
    for i=1:m+1
        if i == 1
            y(x<=pt(i,1)) = pt(i,2);
        elseif i == m+1
            y(x>pt(m,1)) = pt(m,2);
        else
            x1 = pt(i-1,1);
            x2 = pt(i,1);
            y1 = pt(i-1,2);
            y2 = pt(i,2);
            y(x1<x & x<=x2) = (y2-y1)/(x2-x1)*(x(x1<x & x<=x2)-x1)+...
                y1*ones(size(x(x1<x & x<=x2)));
        end
    end
    
    
end
