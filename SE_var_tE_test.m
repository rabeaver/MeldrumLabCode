close all
clc
clear

%%
[ap,spec,spec2] = readTecmag('/Users/tyler/Dropbox/Data/MortarDiffusion/Sample1b_Se_var_tE_nE512_n16_n2D11_24July2014.tnt');
nPts = 35;
oPts = 5;
nEchoes = 512;
oEchoes = 4;
n2D = 11;
tEs = 1e-6*(120:18:300);
tE = zeros(nEchoes-oEchoes,n2D);

for ii = 1:n2D
    tE(:,ii) = (1+oEchoes:nEchoes)*tEs(ii);
end

spec3 = reshape(spec2',nPts,nEchoes,n2D);
spec3 = spec3(1:end-oPts,:,:);
dataInt = abs(sum(real(spec3)));
dataInt = reshape(dataInt,nEchoes,n2D);
dataInt = dataInt(1+oEchoes:end,:);

figure(1)
hold on
for ii = 1:n2D
    dataInt(:,ii) = dataInt(:,ii)./max(dataInt(:,ii));
    plot(tE(:,ii),dataInt(:,ii));
end
hold off

% dataIntI = zeros(nEchoes-oEchoes,n2D);
% 
% for ii = 1:n2D
%     dataIntI(:,ii) = interp1(tE(:,ii),dataInt(:,ii),tE(:,1));
% end

% dataIntI = dataIntI./max(dataIntI,[],1);

% figure(2)
% surf(real(dataIntI)); shading interp
% surf(real(spec3(:,:,1)));

%%
guesses = 0.02;

for ii = 1:n2D
    beta(ii) = nlinfit(tE(:,ii),real(dataInt(:,ii)), @t2monofitonly_simple,guesses); 
end


figure(4)
scatter(tEs, beta)
xlim([0 500e-6]);

cftool(tEs, beta)

