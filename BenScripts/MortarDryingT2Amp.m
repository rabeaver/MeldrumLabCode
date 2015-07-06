clear
close all
clc
% addpath(genpath('Z:\TKM\'));

%%Mortar Drying Matlab
% maindir = 'C:\Users\bmfortman\Documents\Data\';% for the labMACHINE
maindir = 'C:\Users\benjamin\Documents\Data\'; % for my laptop!

dir11 = strcat(maindir,'MortarDrying\OneToOne\1\');
dir21 = strcat(maindir,'MortarDrying\TwoToOne\1\');


%% reading mortar data in then plotting
expNums = [1 : 8];
cd(dir21);
for j = 1:length(expNums)
    
    
    data = strcat('TwoToOne0',num2str(j),'-decaysRe.dat');
    data1(j).sig = load(data);
end
%dataAll = [data1(1).sig; data1(2).sig; data1(3).sig; data1(4).sig; data1(5).sig; data1(6).sig; data1(7).sig; data1(8).sig;];
%expTime = 6.4;
%totTime = expTime*length(dataAll);
guesses = [0.15, 20]%, 0.1, 7];% remember to change guesses for biexp bs monoexp
guesses2 = [0,0]%,0,0];
fitopts = statset('MaxIter',500,'TolX',1e-14,'UseParallel',true,'Display','off');
size = size(data1(1).sig);
% 
% figure(1) %analysis for nick if he needs it
% hold on
% plot(data1(1).sig(:,1),data1(1).sig(:,8))

%% Biexpfit
for i = 1:(size(2)-1)% this is to move along the matrices taking the next point
for j = 1:length(expNums)%this is the exp fit t2bifit needs to be changed to t2monofit to convert between the two
    try [fit(j).beta,fit(j).resid,fit(j).J] = nlinfit((data1(j).sig(:,1)),(data1(j).sig(:,i+1)),@t2monofit_simple,guesses,fitopts);%taking from column 20 for testing
    catch [fit(j).beta,fit(j).resid,fit(j).J] = nlinfit((data1(j).sig(:,1)),(data1(j).sig(:,i+1)),@t2monofit_simple,guesses2,fitopts);%taking from column 20 for testing
    end
    fit(j).pred = t2monofit_simple(fit(j).beta,(data1(j).sig(:,1)));
%     figure(j)
%     hold on % commented out, plotting of data as it is prepared
%     plot((data1(j).sig(:,1)),(data1(j).sig(:,i+1)))%col 20 for testing
%     plot((data1(j).sig(:,1)),fit(j).pred,'-k')
    %[Ypred,delta] = nlpredci(@t2bifit_simple,echoVector',fit(j).beta,fit(j).resid,'Jacobian',fit(j).J);
    ci = nlparci(fit(j).beta,fit(j).resid,'jacobian',fit(j).J);
%     error(j).Ypreds = (sum(Ypred)/(length(Ypred)));
%     error(j).deltas = (sum(delta)/(length(delta)));
    error(j).ci = ci;
    %error(j).delta = delta;
    %bigFit(j).fit = fit.beta;
    %bigError(j).error = error;

end
%big(i).fit = bigFit.fit
%big(i).error = bigError.error
bigFit(i).fit = fit.beta;
bigFit(i).fit2 = fit(2).beta;
bigFit(i).fit3 = fit(3).beta;
bigFit(i).fit4 = fit(4).beta;
bigFit(i).fit5 = fit(5).beta;
bigFit(i).fit6 = fit(6).beta;
bigFit(i).fit7 = fit(7).beta;
bigFit(i).fit8 = fit(8).beta;
bigError(i).error = error;
% bigError(i).error2 = error(2);
% bigError(i).error3 = error(3);
% bigError(i).error4 = error(4);
% bigError(i).error5 = error(5);
% bigError(i).error6 = error(6);
% bigError(i).error7 = error(7);
% bigError(i).error8 = error(8);
% biExp1(i) = fit(1).beta(1);
% biExp1(i) = fit(2).beta(1);
% biExpt2(i) = fit(1).beta(2);
end

%% plotting Biexpfit vs. time and depth
close all
for j = 1:length(bigFit);
    
    bi1Exp1(j) = bigFit(j).fit(1);%The following ones are required for the monoexponential fit
    bi1Expt21(j) = bigFit(j).fit(2);
    bi2Exp1(j) = bigFit(j).fit2(1);
    bi2Expt21(j) = bigFit(j).fit2(2);
    bi3Exp1(j) = bigFit(j).fit3(1);
    bi3Expt21(j) = bigFit(j).fit3(2);
    bi4Exp1(j) = bigFit(j).fit4(1);
    bi4Expt21(j) = bigFit(j).fit4(2);
    bi5Exp1(j) = bigFit(j).fit5(1);
    bi5Expt21(j) = bigFit(j).fit5(2);
    bi6Exp1(j) = bigFit(j).fit6(1);
    bi6Expt21(j) = bigFit(j).fit6(2);
    bi7Exp1(j) = bigFit(j).fit7(1);
    bi7Expt21(j) = bigFit(j).fit7(2);
    bi8Exp1(j) = bigFit(j).fit8(1);
    bi8Expt21(j) = bigFit(j).fit8(2);%end of monoexp fit section
    %beginning of biExp fit section
%     bi1Exp2(j) = bigFit(j).fit(3);
%     bi1Expt22(j) = bigFit(j).fit(4);
%     bi2Exp2(j) = bigFit(j).fit2(3);
%     bi2Expt22(j) = bigFit(j).fit2(4);
%     bi3Exp2(j) = bigFit(j).fit3(3);
%     bi3Expt22(j) = bigFit(j).fit3(4);
%     bi4Exp2(j) = bigFit(j).fit4(3);
%     bi4Expt22(j) = bigFit(j).fit4(4);
%     bi5Exp2(j) = bigFit(j).fit5(3);
%     bi5Expt22(j) = bigFit(j).fit5(4);
%     bi6Exp2(j) = bigFit(j).fit6(3);
%     bi6Expt22(j) = bigFit(j).fit6(4);
%     bi7Exp2(j) = bigFit(j).fit7(3);
%     bi7Expt22(j) = bigFit(j).fit7(4);
%     bi8Exp2(j) = bigFit(j).fit8(3);
%     bi8Expt22(j) = bigFit(j).fit8(4);
    %End of the biexp fit section
%     depth(i).exp1 = biExp.exp1;
%     depth(i).exp2 = biExp.exp2;
%     depth(i).t21 = biExp.t21;
%     depth(i).t22 = biExp.t22;
end
biExp1 = vertcat(bi1Exp1,bi2Exp1,bi3Exp1,bi4Exp1,bi5Exp1,bi6Exp1,bi7Exp1,bi8Exp1);%the concatenation of the monofitExp
biT21 = vertcat(bi1Expt21,bi2Expt21,bi3Expt21,bi4Expt21,bi5Expt21,bi6Expt21,bi7Expt21,bi8Expt21);%concatenation of the T2
% biExp2 = vertcat(bi1Exp2,bi2Exp2,bi3Exp2,bi4Exp2,bi5Exp2,bi6Exp2,bi7Exp2,bi8Exp2);%the concatenation of the biExp
% biT22 = vertcat(bi1Expt22,bi2Expt22,bi3Expt22,bi4Expt22,bi5Expt22,bi6Expt22,bi7Expt22,bi8Expt22);%concatenation of the T2

expDepth = linspace(14500,0,29);
% expTime = 6.4;
% totTime = expTime*1.6704e+3;
% time = linspace(0,totTime,216);
figure(1)
subplot(1,2,1)
hold on
plot(expDepth,biExp1)%mono exp
% plot(expDepth,biExp2,'-r')%biexp
xlabel('Experiment Depth [um]')
ylabel('Amplitude')
legend('1st Biexponent')%, '2nd Biexponent')

subplot(1,2,2)
hold on
axis([0,25000,0,40]);
plot(expDepth,biT21)%monoexp
% plot(expDepth,biT22,'-r')%biexp
xlabel('Experiment Depth[um]')
ylabel('T2 time')
legend('1st T2 time')%, '2nd T2 time')

%% looking at the t2 or amp at a specific depth
% close all

expTime = linspace(0,24,8);

figure(3)
hold on
plot(expTime,biExp1(:,8))
plot(expTime,biExp1(:,12),'-r')
plot(expTime,biExp1(:,16),'-k')
plot(expTime,biExp1(:,24),'-m')
xlabel('time(hours)')
ylabel('amplitude')

figure(4)
hold on
plot(expTime,biT21(:,8))
plot(expTime,biT21(:,12),'-r')
plot(expTime,biT21(:,16),'-k')
plot(expTime,biT21(:,24),'-m')
xlabel('time(hours)')
ylabel('T2 times')
legend('15420[um]','12173[um]','8927[um]','2435[um]')

