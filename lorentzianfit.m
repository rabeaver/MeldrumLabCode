function [y]=fun(v,x)
    gamma = v(1);
    x_0 =   v(2);
    I =     v(3);
    y =     I .* ((gamma.^2)./(4.*(x - x_0).^2 + gamma.^2));
end