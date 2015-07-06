function [y]=fun(v,x)
    I0c =   v(1);
    T2c =   v(2);
    I0i =   v(3);
    T2i =   v(4);
    I0a =   v(5);
    T2a =   v(6);
    nu =    v(7);
    y0 =    v(8);
    
    y =  y0 + I0c.*exp(-(x./T2c).^2).*sin(2*pi*nu.*x)./(2*pi*nu.*x) + I0i.*exp(-x./T2i) + I0a.*exp(-x./T2a);
end


