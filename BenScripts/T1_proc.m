clc
close all
clear all
%% Tecmag Reading in data
cd('C:\Users\benjamin\Documents\MortarDiffusion')%laptop
% cd('C:\Users\bmfortman\Documents\Data\MortarDiffusion')%the nice lab compy
[ap,sp1,sp2] = readTecmag('Jsoakedwater_trep8_T1_n21_nE32_s4_3July2014.tnt');
%% Params and reshaping
nrEchoes = 32;
nrTaupts = 21;
nrPts = 69;
nrBlpts = 5;
nrScans = 4;

data = reshape(sp2',nrPts,32,nrTaupts); % rehapes data into individual echoes
data = data(1:end-nrBlpts,:,:); % cuts out blank points
dataSum = sum(sum(data),2); % sums, each echo, then the echo train

for i = 1:size(dataSum,3) % takes out of a complex 3d matrix, into a linear set of values
    dataSumpts(i) = dataSum(:,:,i);
end
datamax = max(dataSumpts);
dataNorm = dataSumpts./datamax;
linSpacepts = [77519
158993
244839
335554
431723
534043
643357
760690
887314
1024827
1175279
1341360
1526700
1736356
1977687
2261998
2607985
3050039
3662922
4668958
7999994]*1e-6;% this is copy pasted from the t1 estimation excel calculator, can be used for log or linear, just input numbers
%% fitting to t1
guesses = [0.8,.15];
[beta,resid,J] = nlinfit(linSpacepts',real(dataNorm),@T1_recovery_no_y0,guesses);
[pred] = T1_recovery_no_y0(beta,linSpacepts);

figure(1)
hold on
scatter(linSpacepts,real(dataNorm))
plot(linSpacepts,pred)


%% other stuff % these are for troubleshooting the figures to make sure the data is saved correctly
figure(1)
sp2r = real(sp2);
for i=1:21
subplot(3,7,i)
axis([0 2500 0 2e+4])
hold on
plot(sp2r(i,:))
end

figure(2)
surf(sp2r(1:15,:))

figure(3)
surf(sp2r(16:21,:))%
