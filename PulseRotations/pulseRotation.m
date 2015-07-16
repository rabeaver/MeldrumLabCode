function Mt = pulseRotation(M0,t,delw0,ph_d,w1)

                   
ph = ph_d*(2*pi/360);           %rf pulse phase, rad

wnut = sqrt(w1^2+delw0^2);      %nutation frequency

a = 0.5*((w1/wnut)^2+(1+(delw0/wnut)^2)*cos(wnut*t)) + 1i*(delw0/wnut)*sin(wnut*t);
b = 0.5*((w1/wnut)^2+(1+(delw0/wnut)^2)*cos(wnut*t)) - 1i*(delw0/wnut)*sin(wnut*t);
c = (delw0/wnut)^2 + (w1/wnut)^2*cos(wnut*t);
d = (w1/wnut)*((delw0/wnut)*(1-cos(wnut*t)) - 1i*sin(wnut*t))*exp(+1i*ph);
e = (w1/wnut)*((delw0/wnut)*(1-cos(wnut*t)) + 1i*sin(wnut*t))*exp(-1i*ph);
f = (w1/wnut)*((delw0/wnut)*(1-cos(wnut*t)) - 1i*sin(wnut*t))*exp(-1i*ph);
g = (w1/wnut)*((delw0/wnut)*(1-cos(wnut*t)) + 1i*sin(wnut*t))*exp(+1i*ph);
h = 0.5*(w1/wnut)^2*(1-cos(wnut*t))*exp(+1i*2*ph);
j = 0.5*(w1/wnut)^2*(1-cos(wnut*t))*exp(-1i*2*ph);

R = [a h d;
    j b e;
    f g c];

Mt = R*M0;
% Mt = [real(M(1));
%      imag(M(2));
%      real(M(3))];
end