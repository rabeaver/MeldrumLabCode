clear
clc
close all

%%

pos = 0.1:0.05:5;
t = 0.1:0.05:5;

% Adef = zeros(length(pos),1);
% T2def = zeros(length(pos),1);

Adef = exp(-(t-2)/2);
Adef(1:20) = zeros(20,1);

T2def = -0.125*(t-10);
T2def(1:20) = zeros(20,1);

% figure(1)
% hold on
% plot(t,Adef,'-b')
% plot(t,T2def,'-r')

%%

guesses = [0.5 0.5];
fitopts = statset('MaxIter',5000,'TolX',1e-14,'UseParallel',true,'Display','off');

for i = 1:length(pos)
    data(:,i) = Adef(i)*exp(-t/T2def(i)) + 0.05*randn(1,length(pos));
    try
        beta(:,i) = nlinfit(t',data(:,i),@t2monofit_simple,guesses,fitopts);
    catch
        beta(:,i) = [0;0];
    end
    simData(:,i) = t2monofit_simple(beta(:,i),t);
end

posIndex = 30;
pos3D = pos(posIndex)*ones(length(t),1);

% figure(2)
% surf(t,pos,data)
% shading flat
% colormap 'gray'
% xlabel('position')
% ylabel('time')
% zlabel('signal intensity')
% view(210,26)
% plot3(pos3D,t,data(:,posIndex),'-r','LineWidth',2)
% 
% 
% figure(3)
% hold on
% plot(t,t2monofit_simple(beta(:,posIndex),t),'-k')
% plot(t,data(:,posIndex),'-r')
% xlabel('position')
% ylabel('signal intensity')
% ylim([-0.1 1.1])

% figure(4)
% surf(t,pos,simData)
% shading flat
% colormap 'gray'
% xlabel('position')
% ylabel('time')
% zlabel('signal intensity')
% view(210,26)

% figure(5)
% subplot(2,1,1)
% hold on
% plot(pos,beta(1,:),'-r')
% plot(pos,Adef,'-k')
% ylabel('A value (arb)')
% xlabel('position')
% ylim([0 2])
% subplot(2,1,2)
% hold on
% plot(pos,beta(2,:),'-b')
% plot(pos,T2def,'-k')
% ylabel('T_2 value (time)')
% xlabel('position')
% ylim([0 2])


figure(6)
set(6,'Position',[39 508 2429 587])
subplot(2,3,[1 4])
hold on
surf(t,pos,data)
shading flat
colormap 'gray'
xlabel('position')
ylabel('time')
zlabel('signal intensity')
view(210,26)
plot3(pos3D,t,data(:,posIndex),'-r','LineWidth',2)
subplot(2,3,[2 5])
hold on
plot(t,t2monofit_simple(beta(:,posIndex),t),'-k')
plot(t,data(:,posIndex),'-r')
xlabel('time')
ylabel('signal intensity')
ylim([-0.1 1.1])
subplot(2,3,3)
hold on
plot(pos,Adef,'-k')
plot(pos,beta(1,:),'-r')
ylabel('A value (arb)')
xlabel('position')
ylim([0 2])
subplot(2,3,6)
hold on
plot(pos,T2def,'-k')
plot(pos,beta(2,:),'-b')
ylabel('T_2 value (time)')
xlabel('position')
ylim([0 2])

