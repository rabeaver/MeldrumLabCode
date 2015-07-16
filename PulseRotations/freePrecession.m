function Mt = freePrecession(M0,t,delw0)

R = [exp(1i*delw0*t)     0            0;
         0           exp(-1i*delw0*t) 0;
         0               0            1];
 
Mt = R*M0;
end