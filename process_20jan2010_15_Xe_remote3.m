function pd = process_20jan2010_15_Xe_remote3()
%Processes 051008_zslice_xzphase_capillary_velenc_2.fid

%Get raw data with a baseline adjustment based on the last 20 pts of each
% fid.
data=getdata('','20jan2010_15_Xe_remote3.fid',10);

%Reshape the third dimension to reflect the experiment's independent dims.
rawshape=size(data);
data=reshape(data,[rawshape(1:2)]);

%Trim off initial fid pts and apodize.
data = stdfidmanip(data,10,180,[0,0,0,1,1],[0,0,0,10,10],[0,0,0,16.5,16.5],[0,0,0,0,0]);

%FFT the data, first along the direct dimension (can't zero fill yet due to
%memory issues. (Set the echo center to correspond to the original point)
data = fftdata(data,[1,0,0,0,0],[1,1,1,1,1],[rawshape(1:2),3,16,16]);

%Sum over the peak (It's not properly phased but close enough for now) and
%sum over the time of flight data;
data = squeeze(sum(sum(data(100:104,3:17,:,:,:),1),2));

%With a dramatically smaller data set I can now zero fill the image data.
%Must pre zero fill k-space data or not properly center it and get a huge
%phase roll.
data = fftdata(data,[0,1,1],[1,1,1],[3,96,96]);
data = circshift(data,[0,20,60]);

%Calculate the velocity maps based of phase changes
maxmag = max(max(squeeze(abs(data(1,:,:)))));
offset = 1e-6*maxmag;
viewpts = (0.5*maxmag<abs(data(1,:,:)));

basephase = data(1,:,:)./(abs(data(1,:,:))+offset);
data(2,:,:)=viewpts.*data(2,:,:)./basephase;
data(2,:,:)=(180/pi)*atan2(imag(data(2,:,:)),real(data(2,:,:)));
data(3,:,:)=viewpts.*data(3,:,:)./basephase;
data(3,:,:)=(180/pi)*atan2(imag(data(3,:,:)),real(data(3,:,:)));


%Display the data;
figure;
subplot(2,2,1);
imagesc(squeeze(abs(data(1,:,:))));
colorbar;
subplot(2,2,2);
imagesc(squeeze((data(2,:,:))),[0,180]);
colorbar;
subplot(2,2,3);
imagesc(squeeze((data(3,:,:))),[0,180]);
colorbar;

%Send back the results
pd = data;