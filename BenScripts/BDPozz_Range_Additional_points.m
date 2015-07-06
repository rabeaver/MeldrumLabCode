clear
clc
% it clears the data, but not the plots.
%% Doing other fitting

% additional processing, keeping it outside of the script to make it easier
data(2).name = ('Sample0B(1)_CPMG_3262015.tnt');
data(1).name = ('Sample0.5P(1)_CPMG_3262015.tnt');
data(3).name = ('Sample3B(1)_CPMG_3262015.tnt');
data(4).name = ('Sample8B(1)_CPMG_3262015.tnt');
data(5).name = ('Sample15B(1)_CPMG_3262015.tnt');
data(6).name = ('Sample20B(1)_CPMG_3262015.tnt');

dir1 = ('C:\Users\benjamin\Documents\Data\Mortar Curing\B-P Mortar range experiments\'); % personal computer
cd(dir1)

nEchoes = 512
tE = 150e-6
nrPts = 69
nrBlank = 5 

%% fitting and parameters, no plotting in this one
for i = 1:length(data)
    [data(i).params,data(i).data] = readTecmag4d(data(i).name);
    data(i).data = reshape(real(data(i).data),nrPts,nEchoes);
    data(i).data = data(i).data(1:end-nrBlank,:);
    sample(i).data = data(i).data;% conversion, so averaging doesn't have to take place here.
    sample(i).sumData = sum(sample(i).data);
    sample(i).normSum = sample(i).sumData(:)./max(sample(i).sumData);
end

echoAxis = (1:nEchoes)*tE;
guesses = [0.8,0.04,0.2,0.005];
guesses2 = [0.8,0.2];

for i = 1:length(sample)
    [sample(i).beta,sample(i).resid,sample(i).j] = nlinfit(echoAxis(3:end),sample(i).normSum(3:end)',@t2bifit_simple,guesses);
    sample(i).pred = t2bifit_simple(sample(i).beta,echoAxis(3:end));
    sample(i).ci = nlparci(sample(i).beta,sample(i).resid,'jacobian',sample(i).j);
    
    [fixed(i).beta,fixed(i).resid,fixed(i).j] = nlinfit(echoAxis(3:end),sample(i).normSum(3:end)',@t2bifit_monofixed,guesses2);
    fixed(i).pred = t2bifit_monofixed(fixed(i).beta,echoAxis(3:end));
    fixed(i).ci = nlparci(fixed(i).beta,fixed(i).resid,'jacobian',fixed(i).j);

end

%for compiling fits and finding the average
T2 = sample(1).beta(2);
T22 = sample(1).beta(4);
amp = sample(1).beta(1);
amp2 = sample(1).beta(3);
T2eu = sample(1).ci(2,:);
T22eu = sample(1).ci(4,:);
% FT2 = fixed(1).beta(2);
Famp = fixed(1).beta(1);
Fampe = fixed(1).ci(1,:);
Famp2 = fixed(1).beta(2);
Famp2e = fixed(1).ci(2,:);
% FT2eu = fixed(i).ci(2,:);

for i = 2:length(sample)
    T2 = [T2;sample(i).beta(2)];
    T22 = [T22;sample(i).beta(4)];
    amp = [amp; sample(i).beta(1)];
    amp2 = [amp2; sample(i).beta(3)];
    T2eu = [T2eu; sample(i).ci(2,:)];
    T22eu = [T22eu; sample(i).ci(4,:)];
%     FT2 = [FT2;fixed(i).beta(2)];% for fixed values
    Famp = [Famp; fixed(i).beta(1)];
    Fampe = [Fampe; fixed(i).ci(1,:)];
    Famp2 = [Famp2; fixed(i).beta(2)];
    Famp2e = [Famp2e; fixed(i).ci(2,:)];
%     FT2eu = [FT2eu; fixed(i).ci(2,:)];

end

%error stuff to get it in the right format
T2eu = abs([T2eu(:,1) - T2,T2eu(:,2) - T2]);
T22eu = abs([T22eu(:,1) - T22,T22eu(:,2) - T22]);
% FT2eu = abs([FT2eu(:,1) - FT2,FT2eu(:,2) - FT2]);% for fixed fitting 
Fampe = abs([Fampe(:,1) - Famp,Fampe(:,2) - Famp]);
Famp2e = abs([Famp2e(:,1) - Famp2,Famp2e(:,2) - Famp2]);

%% plotting and stuff


figure(1)
scatter(0.5,amp(1),60,'b','x','LineWidth',2)
scatter(0.5,amp2(1),60,'b','x','LineWidth',2)
scatter(0.5,Famp(1),60,'b','x','LineWidth',2)
scatter(0.5,Famp2(1),60,'b','x','LineWidth',2)

scatter([0 3 8 15 20],amp(2:end),60,'r','x','LineWidth',2)
scatter([0 3 8 15 20],amp2(2:end),60,'r','x','LineWidth',2)
scatter([0 3 8 15 20],Famp(2:end),60,'r','x','LineWidth',2)
scatter([0 3 8 15 20],Famp2(2:end),60,'r','x','LineWidth',2)

figure(2)
scatter(0.5,T2(1),60,'b','x','LineWidth',2)
scatter(0.5,T22(1),60,'b','x','LineWidth',2)
scatter([0 3 8 15 20],T2(2:end),60,'r','x','LineWidth',2)
scatter([0 3 8 15 20],T22(2:end),60,'r','x','LineWidth',2)

