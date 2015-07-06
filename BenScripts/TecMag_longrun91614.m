% close all
clear
clc
%% setting up and reading data file and Paramaters
cd('C:\Users\bmfortman\Documents\Data\MortarDrying');
[ap, spec1, spec2] = readTecmag4d('Armory_69pts_256ec_150tE_9-19-14_90-2dpts.tnt');
tEcho = 150e-6;
nEchoes = 256;
nExp = 90;
tRep = 2;
numScans = 256;


echoAxis = (1:nEchoes)*tEcho;% creates echo axis
%% fitting of data
guesses = [0.9,0.03,0.1,0.004];
lowerBounds = [0,0,0,0];
upperBounds = [1,1,1,1];

spec2 = spec2(1:nExp,1:69*nEchoes); 
spec3 = reshape(real(spec2'),69,nEchoes,nExp); % reshaping to more convenient format
data = spec3(1:64,:,:);
dataSum = sum(data); % summing echotrain
dataSumSummed = sum(dataSum,2); % sum of the echotrains

for i = 1:nExp
spec4(i) = dataSumSummed(1,1,i);%converts from 3d to vector
dataNorm(i).echoes = dataSum(1,:,i)./max(dataSum(1,:,i)); % normalizes echoes
[lsqfit(i).beta,lsqfit(i).resid,lsqfit(i).J] =lsqcurvefit(@t2bifit_simple,guesses,echoAxis,dataNorm(i).echoes,lowerBounds,upperBounds);
guesses = lsqfit(i).beta(:);
guesses = guesses';
dataNorm(i).beta = [0,0,0,0];
while dataNorm(i).beta ~= guesses %while loop reallocates guesses
[dataNorm(i).beta,dataNorm(i).R,dataNorm(i).J] = nlinfit(echoAxis,dataNorm(i).echoes,@t2bifit_simple,guesses);
guesses = dataNorm(i).beta;
end
dataNorm(i).pred = t2bifit_simple(dataNorm(i).beta,echoAxis);
t21(i) = dataNorm(i).beta(2);
t22(i) = dataNorm(i).beta(4);
end
spec4 = spec4./max(spec4);

%% Inverse Laplace stuff
alpha = 3e8;% 5e8 for DI water, looking at the other ones 3e8 for Jamestown and Armory both
omitpoints = 1;
% sumData = sum(data1(:,:,i));

lowLim = 1e-3; %min(echoVector)/10000; %
hiLim = 3; %max(echoVector)/10;
nrILTSteps = length(echoAxis)/2; %divides by 2 to make the ilt move faster
for j = 1:nExp
[sample(j).spectrum,sample(j).tau,sample(j).chisq,sample(j).compte] = upnnlsmooth1D(dataNorm(j).echoes',echoAxis',  lowLim, hiLim, alpha ,  -1,  nrILTSteps);
end
iLT = sample(1).spectrum;
for i = 2:nExp
    iLT = [iLT; sample(i).spectrum];% concatenates ilt's
end



%% Presentation of Data
figure(1)
timeAxis = (1:nExp)*numScans * tRep/3600;
plot(timeAxis,spec4)
xlabel('Hours')
ylabel('Relative Amplitude')

figure(3)
hold on
axis([0 14 0 0.04])
plot(timeAxis,t21)
plot(timeAxis,t22,'r')
xlabel('Hours')
ylabel('T_2 times')



figure(2)
hold on
semilogx(sample(1).tau,sample(1).spectrum)