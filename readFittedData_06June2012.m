clear
clc
close all

%%

fid = fopen('/Users/tyler/Dropbox/Painting database/PaintingList2.txt');
C = textscan(fid, '%s %s %u %s %u %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f');
fclose(fid);

Title = C{1};
Artist = C{2};
Year = C{3};
Color = C{4};
ID = C{5};

T1.y0 = C{6};
T1.y0_err = C{7};
T1.A = C{8};
T1.A_err = C{9};
T1.tau = C{10};
T1.tau_err = C{11};
T1.tau_err_p = C{12};

T2.y0 = C{13};
T2.y0_err = C{14};
T2.A = C{15};
T2.A_err = C{16};
T2.tau = C{17};
T2.tau_err = C{18};
T2.tau_err_p = C{19};

T2bi.y0 = C{20};
T2bi.y0_err = C{21};
T2bi.A1 = C{22};
T2bi.A1_err = C{23};
T2bi.tau1 = C{24};
T2bi.tau1_err = C{25};
T2bi.tau1_err_p = C{26};
T2bi.A2 = C{27};
T2bi.A2_err = C{28};
T2bi.tau2 = C{29};
T2bi.tau2_err = C{30};
T2bi.tau2_err_p = C{31};

for i = 1:26
allNumeric(:,i) = C{i+5};
end

T1_numeric = [T1.y0, T1.A, T1.tau];
T2_numeric = [T2.y0, T2.A, T2.tau];
T2_bi_numeric = [T2bi.y0, T2bi.A1, T2bi.tau1, T2bi.A2, T2bi.tau2];

for i  = 1:5
T2_bi_norm(:,i) = T2_bi_numeric(:,i)./std(T2_bi_numeric(:,i));
end

for i = 1:3
    T2_norm(:,i) = T2_numeric(:,i)./std(T2_numeric(:,i));
    T1_norm(:,i) = T1_numeric(:,i)./std(T1_numeric(:,i));
end