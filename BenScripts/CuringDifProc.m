close all 
clear 
clc
%% 
tE = 150e-6;
nEchoes = 512;
nRealpts = 64;
totalPts = 69;

cd('C:\Users\benjamin\Documents\Data\Mortar Curing\CO2 curing experiments\')
[a,postData] = readTecmag4d('ForcecureSample7_11132014_256scans_512Echoes.tnt');
[b,preData] = readTecmag4d('PrecureSample7_1092014_256scans_512Echoes.tnt');

postData = reshape(postData,totalPts,nEchoes);
postData = postData(1:nRealpts,:);
preData = reshape(preData,totalPts,nEchoes);
preData = preData(1:nRealpts,:);

sumPost = sum(postData);
sumPre = sum(preData);

sumPost = real(sumPost./max(sumPost));
sumPre = real(sumPre./max(sumPre));

echoVector = (1:512)*150e-6;

%% Fitting

guesses = [0.8,0.03,0.06,0.005];

[Postbeta,Postresid,Postj] = nlinfit(echoVector,sumPost,@t2bifit_simple,guesses);
[Postpred,Postdelta] = nlpredci(@t2bifit_simple,echoVector,Postbeta,Postresid,'Jacobian',Postj);
PostCi = nlparci(Postbeta,Postresid,'jacobian',Postj);
% PostUci = Postpred+Postdelta;
% PostLci = Postpred-Postdelta;


[Prebeta,Preresid,Prej] = nlinfit(echoVector,sumPre,@t2bifit_simple,guesses);
[Prepred,Predelta] = nlpredci(@t2bifit_simple,echoVector,Prebeta,Preresid,'Jacobian',Prej);
PreCi = nlparci(Prebeta,Preresid,'jacobian',Prej);
% PreUci = Prepred+Predelta;
% PreLci = Prepred-Predelta;

%% creation of matrices from errorbars
Post = [PostCi(:,1),Postbeta',PostCi(:,2)]
Pre = [PreCi(:,1),Prebeta',PreCi(:,2)]

diff = Pre-Post

% formation of vectors for showing difference between the two
%% Plotting

figure(1)
plot(echoVector,Prepred)
hold on
plot(echoVector,Postpred)
legend(strcat('pre',PostStr(1),Post,'Post'))
