clear
clc
close all

% Pulse parameters
delw0 = 10;                      %pulse offset

ph_d = 0;                       %rf pulse phase, degrees 0 = +x, 90 = +y, etc...
ph = ph_d*(2*pi/360);               %rf pulse phase, rad

%pulse length/power cal
flip = 180;                      %flip angle (degrees)
t_cal = 3e-6;                 %pulse length

w1 = sqrt((flip*(2*pi/360)/t_cal)^2-delw0^2); %pulse strength (90)
wnut = sqrt(w1^2+delw0^2);      %nutation frequency

%details for plotting
tp = 3e-6;                       %length of pulse for plotting

tstep = 100;                    %number of time steps

t = linspace(0,tp,tstep);


M_0 = [0;   %x
       0;   %y
       1];  %z

   
% coherence pathway transitions
for i = 1:length(t)
    a = 0.5*((w1/wnut)^2+(1+(delw0/wnut)^2)*cos(wnut*t(i))) + 1i*(delw0/wnut)*sin(wnut*t(i));
    b = 0.5*((w1/wnut)^2+(1+(delw0/wnut)^2)*cos(wnut*t(i))) - 1i*(delw0/wnut)*sin(wnut*t(i));
    c = (delw0/wnut)^2 + (w1/wnut)^2*cos(wnut*t(i));
    d = (w1/wnut)*((delw0/wnut)*(1-cos(wnut*t(i))) - 1i*sin(wnut*t(i)))*exp(+1i*ph);
    e = (w1/wnut)*((delw0/wnut)*(1-cos(wnut*t(i))) + 1i*sin(wnut*t(i)))*exp(-1i*ph);
    f = (w1/wnut)*((delw0/wnut)*(1-cos(wnut*t(i))) - 1i*sin(wnut*t(i)))*exp(-1i*ph);
    g = (w1/wnut)*((delw0/wnut)*(1-cos(wnut*t(i))) + 1i*sin(wnut*t(i)))*exp(+1i*ph);
    h = 0.5*(w1/wnut)^2*(1-cos(wnut*t(i)))*exp(+1i*2*ph);
    j = 0.5*(w1/wnut)^2*(1-cos(wnut*t(i)))*exp(-1i*2*ph);
    
    R(:,:,i) = [a h d;
                j b e;
                f g c];
    
    M_t(:,i) = R(:,:,i)*M_0;
end
    
figure
hold on
plot(t,real(M_t(1,:)),'-k')
plot(t,imag(M_t(2,:)),'-r')
plot(t,real(M_t(3,:)),'-b')
% plot(t,imag(M_t(1,:)),'--k')
% plot(t,imag(M_t(2,:)),'--r')
% plot(t,imag(M_t(3,:)),'--b')
legend('x','y','z')
