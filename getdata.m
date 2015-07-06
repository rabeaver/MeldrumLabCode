function rawdata=getdata(datapath,fidfolder,baselinecorrection)
%Reads the raw data of fidfolder in datapath, modifying it only by a
%baseline correction using the last points in each fid of number
%baselinecorrection.
%The rawdata dimensions increment as... {FID,multi,datablock}
% The user must figure out how to partition the datablock into the
% appropriate dimensions based on the experimental parameters.

%NEED FUNCTIONS
%  readpar, readheader, getblock, 

    %Determine the locations for the various data files.
    nme = fidfolder;
    pth = datapath;
    n = baselinecorrection;
    fullpth = strcat(pth,nme);
    fidname = strcat(fullpth,'/fid');
    parname = strcat(fullpth,'/procpar');
    
    %Determine the dimensioning of the data.
    multi = readpar(parname,'multi'); %number of acquisitions per FID
    at = readpar(parname,'at')/multi; %acquisition time
    
    %Read the data file.
    f = fopen(fidname,'r','b');
    fheader = readheader(f,n);
    ap=fheader.np/2/multi;  %acquired points per FID in the direct dimension
    ap2=fheader.nblocks;    % ???
    
    %--Get first data block
    [d1,bheader] = getblock(f,fheader.np,n);
    d1 = reshape(d1,ap,multi); %Breaks up block in 1st d to read the fid and
                               %  along the 2nd d to read its repeats
                               %  (phase cycles??)
    %--Allocate the overall dataset size
    bshape = size(d1);
    rawdata = zeros([bshape,fheader.nblocks]);
    rawdata(:,:,1)=d1;
    
    %--Get the remaining data blocks
    for fv = 2:(fheader.nblocks)
        [dn,bheader] = getblock(f,fheader.np,n);
        dn = reshape(dn,ap,multi);
        rawdata(:,:,fv)=dn;
    end
    fclose(f);
    
    %Need to figure out how to reshape the rawdata matrix. (They do it
    %manually, make user do it manually latter.
    
end