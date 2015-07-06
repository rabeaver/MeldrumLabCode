% written by Elena; mc.rwth-aachen.de;
% 1. change names of all files !
% 2. check the dimensions (numbers of echoes, number od recovery times...)!
% 3. pay attention to the dimensions if it is 's' or 'ms'.Everything should be in 's'

clear;

R='F:\FLT_Callaghan\test\test_data.txt'; % to read , where it will be red from 
RR='F:\FLT_Callaghan\test\test_'; % to save, where it will be saved to

data_test= load(R);

n=35;                   % number of recovery values
points=5400;              %  number of echoes
correct_cpmg_echoes=points; %correct number of echoes. If the decay is too long and extra 'tail' should be skiped, 
                          % then it is possible to make it shorter (correct_cpmg_echoes < points). 
                          %otherwis keep correct_cpmg_echoes=points.
 
% reading of the data
for i=1:n;           
 for j=1:points;
   Tau_1(i,1)=0.001*data_test(i+(i-1)*points,1); % sec, recovery time
   Tau_2_help(j,1)=0.001*data_test(j,2);  % cpmg, 
                                          % the 0.001 coefficient is to convert from 'ms' to 's',check if you need it
   Data_help(j,i)=data_test(j+(i-1)*points,3);  % amplitude
 end;
end;

%substruction of right echoes
         
for i=1:n;    
 for j=1:correct_cpmg_echoes;         %number of correct cpmg echoes
   Data(j,i)= Data_help(j,i);      
   %Data(j,i)= Data_help(j,i)
   Data_tr_plot(i,j)= Data_help(j,i);     
   Tau_2(j,1)=Tau_2_help(j,1);     
 end;
end;
clear points;
clear n;

[points,n]=size(Data);

% M=Mmax-M look up an article  -> Huerlimann_2002.pdf
for i=1:n;     
  for j=1:correct_cpmg_echoes;
    data_cllghn_I(j,i)=Data(j,n)-Data(j,i);
  end;
end;

for i=1:n-1;           %  minus posledni nulevaja kolonka
 for j=1:correct_cpmg_echoes;     
   Data_minus(j,i)=data_cllghn_I(j,i);
   Data_minus_tr(i,j)=data_cllghn_I(j,i);
   m_Tau1(i,1)=Tau_1(i,1);
 end;
end;

clear TEST_st54;
clear data_test;
clear points;
clear i;
clear j;
clear n;
% the data saving

% just all data
SaveFile=strcat(RR,'data.mat');
save (SaveFile);
clear SaveFile;

% for callaghan
SaveFile=strcat(RR,'data.dat');                      % first file with amplitudes
save ( SaveFile, 'Data_minus','-ASCII');
clear SaveFile;

SaveFile=strcat(RR,'time_1.dat');                    % second file with first time dimension
save ( SaveFile, 'm_Tau1','-ASCII');
clear SaveFile;

SaveFile=strcat(RR,'time_2.dat');                     % third file with second time dimension
save ( SaveFile, 'Tau_2','-ASCII');
clear SaveFile;

plot(Tau_2, Data_help);figure(gcf)
ylabel('signal amplitude [a.u.]','FontSize',12);
xlabel('time [s]','FontSize',12);
savefile=[RR,'cpmg.fig'];
saveas(gcf,savefile);
clear savefile;
savefile=[RR,'cpmg.jpg'];
saveas(gcf,savefile);
clear savefile; 

figure1 = figure('PaperPosition',[0.6345 6.345 20.3 15.23],'PaperSize',[20.98 29.68]);
plot(Tau_1, Data_tr_plot);figure(gcf)
ylabel('signal amplitude [a.u.]','FontSize',12);
xlabel('time [s]','FontSize',12);
savefile=[RR,'sat.fig'];
saveas(gcf,savefile);
clear savefile;
savefile=[RR,'sat.jpg'];
saveas(gcf,savefile);
clear savefile;
