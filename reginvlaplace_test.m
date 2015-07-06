clear
clc
close all

data = load('/Users/tyler/Dropbox/Data/T1T2Testing/T1T2/2/data2D.dat');
T2times = load('/Users/tyler/Dropbox/Data/T1T2Testing/T1T2/2/T2Axis.dat')';

%%

for i = 1:size(data,1)
    [s(:,i),g(:,i),ffit(:,i)] = reginvlaplace(T2times',data(i,:)',size(data,2)-1,0,-1);
end

%%
figure(2)
hold on
for i = 1:size(data,1)
scatter3(log10(s(:,i)),i*ones(size(s,1),1),g(:,i))
end

%%
figure(3)
surf(1:size(data,1),log10(s(:,1)),g)
shading interp