clc
close all
clear 
%% Reads data from text files saved 
% compDir = ('C:\Users\bmfortman\Documents\Data');% lab Compy
% % compDir = ('C:\Users\benjamin\Documents\Data');% personal laptop
% dir1 = strcat(compDir,'\MortarCuring\Brick Dust Experiments');
dir1 = ('C:\Users\bmfortman\Documents\Data\MortarCuring\Brick Dust Experiments')
cd(dir1)
% AX = load('XArm.txt');% These are for Armory the old ones
% AIntY = load('intArm.txt');
% ANormY = load('normArm.txt');
% 
% % DIX = load('XDI.txt');% These are for DI water
% % DInty = load('intDI.txt');
% % DnormY = load('normDI.txt');
% 
% % tX = load('XTap.txt');% for tap water
% % tInty = load('intTap.txt');
% % tnormY = load('normTap.txt');
% 
% JIX = load('XJam.txt');% these are for Jamestown
% JInty = load('intJam.txt');
% JnormY = load('normJam.txt');
% 
% qX= load('Xquik.txt');% for quick set mortar
% qIntY = load('intquik.txt');
% qnormY = load('normquik.txt');
% 
% jwX = load('Xjw.txt');%jamestown water
% jwInt = load('intjw.txt');
% jwNorm = load('normjw.txt');

for i = 1:9
    sample(i).name = strcat('BrickdustSample',num2str(i-1));
% sample(1).name = ('Jamestown');%1-1 samples plus J&A
% sample(2).name = ('Armory');
% sample(3).name = ('Sample1');% nothing added
% sample(4).name = ('Sample2');% Ash only
% sample(5).name = ('Sample3');% brick dust only
% sample(6).name = ('Sample4');% clay only
% sample(7).name = ('Sample5');% all additives
% sample(8).name = ('Sample6');% clay and brick dust only
% sample(9).name = ('Sample13');% 1% brick dust only
% 
% sample(10).name = ('Jamestown');%2-1 samples plus J&A
% sample(11).name = ('Armory');
% sample(12).name = ('Sample7'); %nothing added
% sample(13).name = ('Sample8');% Ash only
% sample(14).name = ('Sample9');% Brick dust only
% sample(15).name = ('Sample10');% Clay only
% sample(16).name = ('Sample11');% all additives
% sample(17).name = ('Sample12');% clay and brick dust only

end
expNums1_1 = (1:9); % for plotting the different figures for each group
% expNums2_1 = (10:17);

expNums = 9; % for all of the experiments
for i = 1:expNums % for reading the data from the text files
    sample(i).xAxis = load(strcat(sample(i).name,'X.txt'));
    sample(i).int = load(strcat(sample(i).name,'int.txt'));
    sample(i).norm = load(strcat(sample(i).name,'norm.txt'));
    sample(i).rho = load(strcat(sample(i).name,'rho.txt'));
    sample(i).t21AV = load(strcat(sample(i).name,'t21AV.txt'));
    sample(i).t22AV = load(strcat(sample(i).name,'t22AV.txt'));
    sample(i).t21A = load(strcat(sample(i).name,'t21A.txt'));
    sample(i).t22A = load(strcat(sample(i).name,'t22A.txt'));
    sample(i).t21x = load(strcat(sample(i).name,'t21x.txt'),'t21x','-ascii');
    sample(i).t22x = load(strcat(sample(i).name,'t22x.txt'),'t22x','-ascii');
    sample(i).t21Ax = load(strcat(sample(i).name,'t21Ax.txt'),'t21Ax','-ascii');
    sample(i).t22Ax = load(strcat(sample(i).name,'t22Ax.txt'),'t22Ax','-ascii');
end
% color = ['b','k','r','c','m','g','y','b','k','r','c','m','g']; % a color vector, so that you can have it change in an iterated loop
color = distinguishable_colors(expNums); % gives a matrix of rgb values for the number of exp's
%%  and then plots it all
figure(1)
title('One to One Distributions')
for i =1% expNums1_1
h1line(i) = semilogx((sample(i).xAxis),sample(i).norm,'color',color(i,:),'LineWidth',2);
axis([0 300 0 1])
hold on

plot(sample(i).t21x,sample(i).t21AV,'color',color(i,:),'LineWidth',2)
plot(sample(i).t22x,sample(i).t22AV,'color',color(i,:),'LineWidth',2)
stem(sample(i).t21Ax,sample(i).t21A,'color',color(i,:))
stem(sample(i).t22Ax,sample(i).t22A,'color',color(i,:))

% semilogx(JIX,JnormY,'-k')
% % semilogx(DIX,DnormY,'-r')
% semilogx(qX,qnormY,'-m')
% % semilogx(tX,tnormY,'-g')
% semilogx(jwX,jwNorm,'--k','LineWidth',2)

end

%% 
xlabel('Pore Diameter [um]')
ylabel('Percent Distribution')
legend(h1line,sample(expNums1_1).name); 
axis([0 120 0 1])
% 
% figure(2)
% title('Two to One Distributions')
% for i =expNums2_1
% h2line(i) = semilogx((sample(i).xAxis),sample(i).norm,'color',color(i,:),'LineWidth',2);
% axis([0 300 0 1])
% hold on
% 
% plot(sample(i).t21x,sample(i).t21AV,'color',color(i,:),'LineWidth',2)
% plot(sample(i).t22x,sample(i).t22AV,'color',color(i,:),'LineWidth',2)
% stem(sample(i).t21Ax,sample(i).t21A,'color',color(i,:))
% stem(sample(i).t22Ax,sample(i).t22A,'color',color(i,:))
% 
% % semilogx(JIX,JnormY,'-k')
% % % semilogx(DIX,DnormY,'-r')
% % semilogx(qX,qnormY,'-m')
% % % semilogx(tX,tnormY,'-g')
% % semilogx(jwX,jwNorm,'--k','LineWidth',2)
% 
% end
% xlabel('Pore Diameter [um]')
% ylabel('Percent Distribution')
% % legend(h2line,sample(expNums2_1).name); 
% 
% % legend(h2line,sample(1:expNums).name); 
% axis([0 120 0 1])



figure(3)
title('Integrated Data')
for i = expNums1_1;
semilogx((sample(i).xAxis),sample(i).int,'color',color(i,:))
axis([0 300 0 1])
hold on
end
% semilogx(JIX,JInty,'-k')
% % semilogx(DIX,DInty,'-r')
% semilogx(qX,qIntY,'-m')
% semilogx(jwX,jwInt,'-.k')
xlabel('Pore Diameter [um]')
ylabel('Percent Distribution')
% legend(sample(1:expNums).name);
axis([0 120 0 1])
% 
% figure(4)
% title('Integrated Data')
% for i = expNums2_1;
% semilogx((sample(i).xAxis),sample(i).int,'color',color(i,:))
% axis([0 300 0 1])
% hold on
% end
% % semilogx(JIX,JInty,'-k')
% % % semilogx(DIX,DInty,'-r')
% % semilogx(qX,qIntY,'-m')
% % semilogx(jwX,jwInt,'-.k')
% xlabel('Pore Diameter [um]')
% ylabel('Percent Distribution')
% % legend(sample(1:expNums).name);
% axis([0 120 0 1])
% 
% axis([1 100 0 1]); xlabel('Diameter [um]'); ylabel('Percent Abundance');

%% surface plot (pore size)
spectrum = sample(1).norm;
for i = 2:9
    spectrum = [spectrum;sample(i).norm]
end
figure(4)

surf((sample(1).xAxis),(expNums1_1-1)*0.25,spectrum)
xlabel('Pore Diameter [um]')
xlim([0 40])
ylabel('Brick Dust %')
zlabel('Percent Distribution')
title('Samples [1-8] with increasing brick Dust from 0.25 - 2 percent')
shading flat