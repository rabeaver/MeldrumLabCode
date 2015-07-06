function [y]=fun(v,x)
    gamma1 = v(1);
    x_01 =   v(2);
    I1 =     v(3);
    gamma2 = v(4);
    x_02 =   v(5);
    I2 =     v(6);
    gamma3 = v(7);
    x_03 =   v(8);
    I3 =     v(9);
    
    y =     I1 .* ((gamma1.^2)./(4.*(x - x_01).^2 + gamma1.^2)) + I2 .* ((gamma2.^2)./(4.*(x - x_02).^2 + gamma2.^2)) + I3 .* ((gamma3.^2)./(4.*(x - x_03).^2 + gamma3.^2));
end