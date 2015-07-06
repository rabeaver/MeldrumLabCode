clear
clc
close all
% 
% data = load('data2D.dat');

fileid = fopen('GlyerolWaterGd.pt2');

% owner = fread(fileid,4,'int8=>char'); 
% format = fread(fileid,4,'int8=>char'); 
% nrSubplotRows = fread(fileid,1,'int16'); 
% nrSubplotCols = fread(fileid,1,'int16'); 
% version = fread(fileid,4,'int8=>char'); 
% 
% % fseek(fileid,874,'cof');
% 
% xMajTickLength = fread(fileid,1,'float'); 
% xMinTickLength = fread(fileid,1,'float'); 
% xTickSpac= fread(fileid,1,'float'); 
% yTickSpac= fread(fileid,1,'float'); 
% xTicksPerLabel= fread(fileid,1,'float'); 
% yTicksPerLabel= fread(fileid,1,'float'); 
% 
% axesFont = fread(fileid,60,'*char');
% labelFont = fread(fileid,60,'*char');
% titleFont = fread(fileid,60,'*char');
% 
% axesFontColor = fread(fileid,4,'*char');
% labelFontColor = fread(fileid,4,'*char');
% titleFontColor = fread(fileid,4,'*char');
% 
% yLabelVert = fread(fileid,1,'bool');
% axesMode = fread(fileid,1,'short');
% drawXGrid = fread(fileid,1,'bool');
% drawYGrid = fread(fileid,1,'bool');
% drawXFine = fread(fileid,1,'bool');
% drawYFine = fread(fileid,1,'bool');
% 
% plotTitle = fread(fileid,1,'*char');
% 
% colorMapLength = fread(fileid,4,'int32=>long'); 
% pos = ftell(fileid);
fseek(fileid,-106496,'eof');
data = fread(fileid,106496,'float');

% xdim = fread(fileid,4,'int32=>long');
% 
% dataType = fread(fileid,1,'int32'); 
% xDim = fread(fileid,1,'int32'); 
% yDim = fread(fileid,1,'int32'); 
% zDim = fread(fileid,1,'int32'); 
% qDim = fread(fileid,1,'int32'); 




fclose(fileid);

data = reshape(data,2048,13);

data;