function [processedSpectrum, processedFID, ph0, regFID] = processFID(fidfolder,datapath,lsfid,baselinepts,lb,multi,sw,at,phasenumber,bcstart,bcstop)
    
    nme = fidfolder;
    pth = datapath;
    n = baselinepts; %baselinecorrection; %use this for remote acquisitions with "multi"
    fullpth = strcat(pth,nme);
    fidfilename = strcat(fullpth,'/fid');
    parfilename = strcat(fullpth,'/procpar');

    %Read the data file.
    
    f = fopen(fidfilename,'r','b');
    fheader = readheader(f,n);
    ap=fheader.np/2;  %acquired points per FID in the direct dimension -- complex points
    ap2=fheader.nblocks;    % number of FIDs

    
    %--Get first data block
    [d1,bheader] = getblock(f,fheader.np,n);
    
    
    %use this for remote acquisitions with "multi"
    d1 = reshape(d1,ap/multi,multi); %Breaks up block in 1st d to read the fid and along the 2nd d to read its repeats (phase cycles??)
    
    %--Allocate the overall dataset size
%     bshape = size(d1,1);
    bshape = size(d1); %use this for remote acquisitions with "multi"
%     rawdata = zeros([bshape,fheader.nblocks]);
    rawdata(:,:,1)=d1; %use this for remote acquisitions with "multi"
%     rawdata(:,1)=d1;
    

    
    %--Get the remaining data blocks
    for fv = 2:(fheader.nblocks)
        [dn,bheader] = getblock(f,fheader.np,n);
        dn = reshape(dn,ap/multi,multi); %use this for remote acquisitions with "multi"
        rawdata(:,:,fv)=dn; %use this for remote acquisitions with "multi"
        %rawdata(:,fv)=dn
    end
    fclose(f);
    
    rawdata2 = reshape(rawdata,ap,ap2); %change ap2 to "1" for 1D TOF curves
    regFID = rawdata2;
 
    %%%%%%%%%%

    % move on to processing the spectra

    varianFID = rawdata;
    numberFIDs = size(varianFID,2)*fheader.nblocks;

    sizeFID = size(varianFID,1);
    fidpoints = (0:1:sizeFID-1)';

    %weighting for line broadening
    weighting = exp(-0.5.*fidpoints*pi*(lb/sw));

    % left shift data
    varianFID = varianFID(lsfid+1:sizeFID,:);
    varianFID(sizeFID-lsfid:sizeFID,:) = 0;
    
    % baseline correct
    realFID = real(varianFID);
    imagFID = imag(varianFID);

    realFIDbc = realFID(sizeFID - baselinepts:sizeFID,:);
    imagFIDbc = imagFID(sizeFID - baselinepts:sizeFID,:);

    meanReal = mean(realFIDbc);
    meanImag = mean(imagFIDbc);

    
    for i=1:1:numberFIDs
        %baseline correction
        realFID(:,i) = realFID(:,i) - meanReal(i);
        imagFID(:,i) = imagFID(:,i) - meanImag(i);
        
        % apodize to a cosine and apply linebroadening
        realFID_sin(:,i) = realFID(:,i).*(cos((pi/2)*fidpoints./sizeFID)).*weighting;
        imagFID_sin(:,i) = imagFID(:,i).*(cos((pi/2)*fidpoints./sizeFID)).*weighting;

        %recombine real and imag parts
        varianFID(:,i) = (realFID_sin(:,i) + 1i*imagFID_sin(:,i));
    end

    % ZERO-FILL
    varianFID(sizeFID+1:2^(nextpow2(varianFID)+1),:) = 0;
    processedFID = varianFID;

    % FT
    varianSpectrum = fftshift(fft(varianFID),1);

    % phasing
    % auto phase to most first signal
    [C,I] = max(abs(varianSpectrum));
    phasedSpectrum = varianSpectrum .* exp(-1i*angle(varianSpectrum(I(phasenumber))));
    ph0 = 180*angle(varianSpectrum(I(phasenumber)))/pi;
    
    % baseline correct spectra
    realSpectrum = real(phasedSpectrum);
    imagSpectrum = imag(phasedSpectrum);

    realSpectrumbc = realSpectrum(bcstart:bcstop,:);
    imagSpectrumbc = imagSpectrum(bcstart:bcstop,:);

    meanReal = mean(realSpectrumbc);
    meanImag = mean(imagSpectrumbc);

    
    for i=1:1:numberFIDs
        %baseline correction
        realSpectrum(:,i) = realSpectrum(:,i) - meanReal(i);
        imagSpectrum(:,i) = imagSpectrum(:,i) - meanImag(i);
        
        %recombine real and imag parts
        phasedSpectrumBC(:,i) = (realSpectrum(:,i) + 1i*imagSpectrum(:,i));
    end
    
    processedSpectrum = phasedSpectrumBC;

end


