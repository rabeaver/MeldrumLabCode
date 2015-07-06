%Get raw data with a baseline adjustment based on the last 20 pts of each
% fid.
%data=getdata('/Users/tyler/Documents/Research/Data/RemoteApr2010/','13apr2010_1_Xe_travelcurve.fid',10);
data=getdata('/Users/tyler/Documents/Research/Data/RemoteApr2010/','01apr2010_1_Xe_travelcurve.fid',10);

%Reshape the third dimension to reflect the experiment's independent dims.
rawshape=size(data);
data=reshape(data,[rawshape(1:2)]);

%Trim off initial fid pts and apodize.
data = stdfidmanip(data,10,180,[0,0,0,1,1],[0,0,0,10,10],[0,0,0,16.5,16.5],[0,0,0,0,0]);

%FFT the data, first along the direct dimension (can't zero fill yet due to
%memory issues. (Set the echo center to correspond to the original point)
data = fftdata(data,[1,0,0,0,0],[1,1,1,1,1],[rawshape(1:2),3,16,16]);

%Display the data;
multi = 60;
plotdim = round(sqrt(multi));

figure
% plotindex1 = 1;
% plotindex2 = 1;
for j=1:1:multi
    subplot(plotdim,plotdim,j)
    plot(abs(data(:,j)));
%     plotindex2 = plotindex2 + 1;
%     if plotindex2 == plotdim1
%         plotindex2 = 1;
%         plotindex1 = plotindex1 + 1;
%     end
end
