clear
close all
clc
%% Start of data analysis for mortar masses
% started on 10/16/14 moved on
masses = [22.8843 22.8942 22.8842 22.8795 22.8725; %Armory sample
    55.5952 55.6139 55.6034 55.6044 55.5784;% Jamestown Sample
    46.2263 46.2248 46.2252 46.2123 46.1923;%Sample 112 
    45.0362 45.0326 44.9946 44.9762 44.9559;% Sample 52
    46.8726 46.8726 46.8613 46.8415 46.8233]; % Sample 721

time = [1:length(masses)]; % computes a time vector

massesDif = diff(masses,1,2);% computes the differences element by element
%% plotting the figures
color = ['b','k','r','g','c'];
figure(1)
for i = 1:size(masses,2)
    plot(time,masses(i,:),color(i))
    hold on
end
legend('Armory',' Jamestown','Sample 112','Sample52','Sample721')
title('masses')

figure(2)
for i = 1:size(masses,2)
    plot(time(2:end),massesDif(i,:),color(i))
    hold on
end
legend('Armory',' Jamestown','Sample 112','Sample52','Sample721')
title('rate of change of masses')
