%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                        %
%   Read binary file writtten in Yi-Qiao Song's format.                                  %
%   The first 20 lines contain the header information. Most of this is disregarded       %
%   except for the echo spacing, the number of recovery times and the number of          %
%   echoes in the CPMG.                                                                  %
%                                                                                        %
%   The function reads the data and recovery times from the specified files and writes   %
%   the data, tau_1 and tau_2 and NoiseStd values to specified file.                     %
%                                                                                        %
%   The Noise std is computed as the standard deviation of the noise in the              %
%   imaginary channel (after phase rotation). The data are stored as complex numbers     %
%   in Yi-Qiao's format.                                                                 %
%                                                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


DataDir = 'crudeoils';
DataFileName = 'X2.ser.ar'; 
Tau_1FileName = 'vdlist'; 

WriteFileName = strcat(DataFileName,'.mat');

cd(DataDir);

% Open the file containing the Tau_1 values
fid = fopen(Tau_1FileName, 'r');
if (fid ~= -1) [Tau_1, count] = fscanf(fid, '%g', inf);
else fprintf(1, 'Error opening file %s. Does it exist ??\n', Tau_1FileName); return; end
fclose(fid);





%Open the file containing the data
fid = fopen(DataFileName,'r');
if (fid == -1) fprintf(1, 'Error opening file %s. Does it exist ??\n', DataFileName); return; end;

%Read the header and get info about echo-spacing, number of wait times and number of echoes
i = 1;
%EndOfBuffer = 0;
%while (~EndOfBuffer)
EndOfBuffer = [];
while( isempty(EndOfBuffer))
   line = fgetl(fid);
   fprintf('%d line = %s\n', i,line);
   
   % Find number of recovery times
   k = findstr(line, 'al1');
   k1 = isempty(k);
   k2 = not(k1);
   if (k2)
      k1 = findstr(line, '='); 
      N_Tau1 = str2num(line(k1+1:end)); 
   end
   
   % Find number of echoes in CPMG sequence
   k = findstr(line,'al2');
   if (~isempty(k)) k1 = findstr(line, '='); N_Tau2 = str2num(line(k1+1:end)); end
   
   % Find the echo spacing 
   k = findstr(line, 'dwell2');
   if (~isempty(k))  
      k1 = findstr(line, '='); 
      k2 = findstr(line, 's');
      Te = str2num(line(k1+1:k2-1));
   end
   
   i = i+1;
   %EndOfBuffer = strcmp(line, 'end_of_header');
   EndOfBuffer = findstr(line, 'end_of_header');
end

fprintf(1, 'Reading header\n');
%Read the data until the end of the file
fprintf(1, 'Reading data ... \n');
[A, count] = fread(fid, inf,'float32'); fclose(fid);

% Every second number belongs to the imaginary channel
Noise = A(2:2:end);
% Alternate numbers from the beginning belong to the real channel
A = A(1:2:end);
%Compute the number of wait times
Number_Tau_1 = length(Tau_1);
% Compute the number of echoes
Number_Tau_2 = length(A)/Number_Tau_1;

if ( (Number_Tau_1 ~= N_Tau1) | (Number_Tau_2 ~= N_Tau2))
   fprintf(1, 'Error in number of Tau1s \n');
   fprintf(1, 'Make sure that you have the right set of files\n');
else
   Tau_2 = [Te : Te : Number_Tau_2 * Te]';
   %Write data from array A into data matrix.
   clear Data
   for i  = 1:Number_Tau_1
      Data(i,:) = A((i-1)*Number_Tau_2 + 1 : i*Number_Tau_2)';
   end
   
   Data = Data';
   NoiseStd = std(Noise);
   
   figure
   mesh(Tau_1, Tau_2, Data);
   xlabel('Tau_1')
   ylabel('Tau_2')
   zlabel('Data')
   save(WriteFileName,'Data' ,'Tau_1', 'Tau_2', 'NoiseStd');
   fprintf(1, 'Data saved in file %s \n', WriteFileName);
end
