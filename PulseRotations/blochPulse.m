clear
close all
clc

%%

M0 = [0;    %Mx
      0;    %My
      1];   %Mz
  
a = pi/2;   %flip angle

% Rotation matrices
Rz = [ cos(a)  -sin(a)  0      ;
       sin(a)  cos(a)   0      ;
       0       0        1      ];
   
Rx = [ 1       0        0      ;
       0       cos(a)   sin(a) ;
       0       -sin(a)  cos(a) ];
   
Ry = [ cos(a)  0        sin(a) ;
       0       1        0      ;
       -sin(a) 0        cos(a) ];
   
M = Ry*M0   

    