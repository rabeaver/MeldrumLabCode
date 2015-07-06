close all
clc
clear

%%
dir = 'C:\Users\tkmeldrum\Dropbox\Data\MortarDrying\OneToOne\1\';
expNum = 8;
fileStem = 'OneToOne';
cd(dir)

for i=1:expNum
    data(i).dat = load(strcat(fileStem,'0',num2str(i),'.dat'));
    data(i).Re = load(strcat(fileStem,'0',num2str(i),'-decaysRe.dat'));
    data(i).Im = load(strcat(fileStem,'0',num2str(i),'-decaysIm.dat'));
    data(i).Cp = complex(data(i).Re,data(i).Im);
end

expTime = load(strcat(fileStem,'exptime.dat'));


figure(1)
hold on
for i = 1:expNum
    plot(data(i).dat(:,1),data(i).dat(:,2))
end