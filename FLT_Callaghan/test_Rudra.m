clear;

for i=1:64;
    Tau_1(i,1)=5*i; 
    Tau_2(i,1)=1*i;
end;

for i=1:64;
for j=1:64;
Data(i,j)=(exp(-Tau_1(i,1)/100))*exp(-Tau_1(j,1)/50);
end;
end;

surf (Data)
NoiseStd=120;

clear i j

 save L:\Aachen_RWTH\FLT\test.mat Data Tau_1 Tau_2 NoiseStd
  save L:\Aachen_RWTH\FLT\test_M Data -ASCII;
  save L:\Aachen_RWTH\FLT\test_t1 Tau_1 -ASCII;
  save L:\Aachen_RWTH\FLT\test_t2 Tau_2 -ASCII;

 FlI