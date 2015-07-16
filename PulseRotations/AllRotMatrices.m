M0 = [0; 0; 1];

a = pi/2;

Rx = [1 0 0;
      0 cos(a) -sin(a);
      0 sin(a) cos(a)];
Ry = [cos(a) 0 sin(a);
      0     1   0;
      -sin(a) 0 cos(a)];
  
Rz = [cos(a) -sin(a) 0;
      sin(a) cos(a)  0;
      0     0       1];

Mx = Rx*M0;
My = Ry*Mx;
Mz = Rz*My