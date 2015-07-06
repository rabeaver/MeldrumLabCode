function [ap,spec,spec2,spec3] = readTecmag(pathname)
% -----------------------------------------------------------------
% Read in Tec Mag NTNMR data
% -----------------------------------------------------------------
% DATATYPE            = definedatatype;
% ap                  = initacquparm;
% ap.dtype            = DATATYPE.TecMag;      % TecMag NTNMR data
spec                = 0;
spec2               = 0;
spec3               = 0;
% -----------------------------------------------------------------
% Read the NTNMR file header 
% -----------------------------------------------------------------
[fileid,message]    = fopen(pathname,'r','ieee-le');
if fileid == -1 
    warndlg(message,'Open Data File');
    return
end

paLen               = 20;                           % preamble length
fseek(fileid,16,'bof');
tsLen               = fread(fileid,1,'int32');      % tsLen of TecMag structure

% --- start of Tecmac struct
td1                 = fread(fileid,1,'int32');      % requested points, npts[1] (t2, acquisition)
td2                 = fread(fileid,1,'int32');      % npts[2] (t1)
td3                 = fread(fileid,1,'int32');      % npts[3]
td4                 = fread(fileid,1,'int32');      % npts[4]
atd1                 = fread(fileid,1,'int32');      % requested points, npts[1] (t2, acquisition)
atd2                 = fread(fileid,1,'int32');      % npts[2] (t1)
atd3                 = fread(fileid,1,'int32');      % npts[3]
atd4                 = fread(fileid,1,'int32');      % npts[4]
fseek(fileid,paLen+52,'bof');
ns                  = fread(fileid,1,'int32');      % number of scans 1D requested, scans
fseek(fileid,paLen+84,'bof');
sfo1                = fread(fileid,1,'float64');    % spectrometer frequency, ob_freq[1] 
sfo2                = fread(fileid,1,'float64');    % spectrometer frequency, ob_freq[2]
sfo3                = fread(fileid,1,'float64');    % spectrometer frequency, ob_freq[3]
sfo4                = fread(fileid,1,'float64');    % spectrometer frequency, ob_freq[4]
fseek(fileid,paLen+180,'bof');
ref                 = fread(fileid,1,'float64');    % reference frequency, ref_freq
                      fread(fileid,1,'float64');    % NMR_frequency
obs_channel         = fread(fileid,1,'int16');      % observe channel
fseek(fileid,paLen+240,'bof');
swh1                = fread(fileid,1,'float64');    % spectral width in Hz, sw[1] (t2, acquisition)
swh2                = fread(fileid,1,'float64');    % sw[2] (t1)
swh3                = fread(fileid,1,'float64');    % sw[3]
swh4                = fread(fileid,1,'float64');    % sw[4]
dw1                 = fread(fileid,1,'float64');    % dwel1[1]
dw2                 = fread(fileid,1,'float64');    % dwel1[2]
dw3                 = fread(fileid,1,'float64');    % dwel1[3]
dw4                 = fread(fileid,1,'float64');    % dwel1[4]
fseek(fileid,paLen+864,'bof');
date                = fread(fileid,32,'*char');     % date[]
nuc1                = fread(fileid,16,'*char');     % nucleus[]
nuc2                = fread(fileid,16,'*char');     % nucleus_2D[]
nuc3                = fread(fileid,16,'*char');     % nucleus_3D[]
nuc4                = fread(fileid,16,'*char');     % nucleus_4D[]
% --- end of Tecmac struct

ctd                 = td1;
date                = strtok(date,0);
nuc1                = strtok(nuc1,0);
nuc2                = strtok(nuc2,0);
nuc3                = strtok(nuc3,0);
nuc4                = strtok(nuc4,0);

swh1                = swh1 * 2;                 % spectral width in hz, sw +/- on TecMag
swh2                = swh2 * 2;                 % spectral width in hz, sw +/- on TecMag
swh3                = swh3 * 2;                 % spectral width in hz, sw +/- on TecMag
swh4                = swh4 * 2;                 % spectral width in hz, sw +/- on TecMag

% -----------------------------------------------------------------
% Set up the acquisition parameters
% -----------------------------------------------------------------
ap.ctd              = ctd;                      % complex points in the time domain
ap.date             = date';                    % date
ap.dw(1)            = dw1;                      % dwell time (interval between pts in the time domain)
ap.dw(2)            = dw2;                      % dwell time (interval between pts in the time domain)
ap.dw(3)            = dw3;                      % dwell time (interval between pts in the time domain)
ap.dw(4)            = dw4;                      % dwell time (interval between pts in the time domain)
ap.ns               = ns;                       % number of scans
ap.nuc1             = deblank(nuc1');           % nucleus #1
ap.nuc2             = deblank(nuc2');           % nucleus #2
ap.nuc3             = deblank(nuc3');           % nucleus #3
ap.nuc4             = deblank(nuc4');           % nucleus #4
ap.pulprog          = 'ntnmr';                  % pulseprogram
ap.swh(1)           = swh1;                     % spectral width in Hz
ap.swh(2)           = swh2;                     % spectral width in Hz
ap.swh(3)           = swh3;                     % spectral width in Hz
ap.swh(4)           = swh4;                     % spectral width in Hz
ap.td(1)            = td1;                      % points in the time domain 
ap.td(2)            = td2;                      % points in the second time domain
ap.td(3)            = td3;                      % points in the third time domain
ap.td(4)            = td4;                      % points in the fourth time domain
ap.aq               = ap.dw(1) * (ctd-1);       % acquisition time (approximate)
ap.atd(1)           = atd1;
ap.atd(2)           = atd2;
ap.atd(3)           = atd3;
ap.atd(4)           = atd4;

% -----------------------------------------------------------------
%  parmode: 1 = 1D, 2 = 2D, 3 = 3D
% -----------------------------------------------------------------
parmode             = 1;
if td2 > 1
    parmode = parmode + 1;
    if td3 > 1
        parmode = parmode + 1;
        if td4 > 1
            parmode = parmode + 1;
        end
    end
end
ap.parmode          = parmode;

% -----------------------------------------------------------------
%  the frequency of the obs_channel needs to be associated 
%  with the first experimental dimension
% -----------------------------------------------------------------
if obs_channel == 3
    ap.sfo(1)   = sfo3;
    ap.sfo(2)   = sfo2;
    ap.sfo(3)   = sfo1;
elseif obs_channel == 2
    ap.sfo(1)   = sfo2;
    ap.sfo(2)   = sfo1;
    ap.sfo(3)   = sfo3;
else
    ap.sfo(1)   = sfo1;                         % spectrometer frequency # 1
    ap.sfo(2)   = sfo2;                         % spectrometer frequency # 2
    ap.sfo(3)   = sfo3;                         % spectrometer frequency # 3
end
ap.sfo(4)           = sfo4;                     % spectrometer frequency # 4

% -----------------------------------------------------------------
%  calculated parameters
% -----------------------------------------------------------------
if ap.sfo(1) > 0.0
    ap.sw(1)        = ap.swh(1) / ap.sfo(1);    % sweep width in ppm
    ap.offset(1)    = ap.sw(1) / 2 + (ref / ap.sfo(1)); % left edge of spectrum (for plotting)
else
    ap.sw(1)        = ap.swh(1);
    ap.offset(1)    = ap.sw(1) / 2;             % left edge of spectrum (for plotting)
end
ap.reshz(1)         = ap.swh(1) / (ctd-1);      % digital resolution in Hz
ap.resppm(1)        = ap.sw(1)  / (ctd-1);      % digital resolution in ppm
if td2 > 1
    if ap.sfo(2) > 0.0                          % hetero nuclear (?)
        ap.sw(2)    = ap.swh(2) / ap.sfo(2);    % sweep width in ppm
    else                                        % homonuclear (?)
        %ap.swh(2)  = ap.swh(1);                % sweep width in Hz
        ap.sfo(2)   = ap.sfo(1);                % spectrometer frequency # 2
        ap.sw(2)    = ap.sw(1);                 % sweep width in ppm
        ap.nuc2     = ap.nuc1;                  % nucleus
    end
    ap.offset(2)    = ap.sw(2) / 2;             % left edge of spectrum (for plotting)
    ap.reshz(2)     = ap.swh(2) / (td2-1);      % digital resolution in Hz
    ap.resppm(2)    = ap.sw(2)  / (td2-1);      % digital resolution in ppm
end
if td3 > 1
    if ap.sfo(3) > 0.0
        ap.sw(3)    = ap.swh(3) / ap.sfo(3);
    else
        ap.sw(3)    = ap.swh(3);
    end
    ap.offset(3)    = ap.sw(3) / 2;             % left edge of spectrum (for plotting)
    ap.reshz(3)     = ap.swh(3) / (td3-1);      % digital resolution in Hz
    ap.resppm(3)    = ap.sw(3)  / (td3-1);      % digital resolution in ppm
end
if td4 > 1
    if ap.sfo(4) > 0.0
        ap.sw(4)    = ap.swh(4) / ap.sfo(4);
    else
        ap.sw(4)    = ap.swh(4);
    end
    ap.offset(4)    = ap.sw(4) / 2;             % left edge of spectrum (for plotting)
    ap.reshz(4)     = ap.swh(4) / (td4-1);      % digital resolution in Hz
    ap.resppm(4)    = ap.sw(4)  / (td4-1);      % digital resolution in ppm
end

% -----------------------------------------------------------------
% Read the NTNMR file data
% -----------------------------------------------------------------
fseek(fileid,paLen+tsLen+12,'bof');
if td2 == 1 && td3 == 1
    [a,cnt] = fread(fileid,td1*2,'float32');
    spec    = a(1:2:cnt) - j*a(2:2:cnt);
    if cnt < td1
        spec(cnt/2:ctd)=0;
    end
    spec2   = 0;
elseif td2 > 1 && td3 == 1
    spec2   = zeros(td2,td1);
    for k=1:td2
        [a,cnt]     = fread(fileid,[1 td1*2],'float32');
        % spec2(k,:)  = a(1:2:cnt) - j*a(2:2:cnt);
        if cnt == td1*2
            spec2(k,:)  = a(1:2:cnt) - j*a(2:2:cnt);
%         else
%             k
%             break;
        end
    end
    spec    = spec2(1,:).';
elseif td2>1 && td3>1
    spec3 = zeros(td3,td2,td1);
    for l =1:td3
        for k=1:td2
            [a,cnt]     = fread(fileid,[1 td1*2],'float32');
            % spec2(k,:)  = a(1:2:cnt) - j*a(2:2:cnt);
            if cnt == td1*2
                spec3(l,k,:)  = a(1:2:cnt) - j*a(2:2:cnt);
                %         else
                %             k
                %             break;
            end
        end
    end
end

% --- start of Tecmac2 struct
daLen       = 2*4*td1*td2;
fseek(fileid,paLen+tsLen+12+daLen+12+828,'bof');
ap.smode    = fread(fileid,1,'int16');     % spectrum mode (0 - tdom; 1 fdom)
% --- end of Tecmac2 struct

fclose(fileid);

