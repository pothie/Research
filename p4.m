function y = p4(x)
    y = zeros(size(x));
    a = 0.5;
    z = -0.7;
    d = 0.005;
    l = 10;
    b = log(2)/(36*d^2);
    G =@(x,b,z) exp(-b.*(x-z).^2);
    F =@(x,l,a) sqrt(max([1-l^2*(x-a).^2;zeros(1,length(x))]));
    A = -0.8<=x&x<=-0.6;
    y(A)=(1/6)*(G(x(A),b,z-d)+G(x(A),b,z+d)+4*G(x(A),b,z));
    B = -0.4<=x&x<=-0.2;
    y(B) = 1;
    C = 0<=x&x<=0.2;
    y(C) = 1-abs(10*(x(C)-0.1));
    D = 0.4<=x&x<=0.6;
    y(D)=(1/6)*(F(x(D),l,a-d)+F(x(D),l,a+d)+4*F(x(D),l,a));
 end