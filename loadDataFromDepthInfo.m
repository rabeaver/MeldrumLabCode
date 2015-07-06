function [params,measDepth,depth,T1time,allT1Re,allT1Im,allT1Mag,T2time,allT2Re,allT2Im,allT1Phased,allT2Phased,T2rough,T2roughPhased,T2phaseAngle] = loadDataFromDepthInfo(dataDir,parfilestem,autoPhase)
    startDir = cd;
    
    cd(dataDir)
    
    %%%
    params.acqTime = readpar_Kea(strcat(parfilestem,'.par'),'acqTime');
    params.bandwidth = readpar_Kea(strcat(parfilestem,'.par'),'bandwidth');
    params.nrScansT1 = readpar_Kea(strcat(parfilestem,'.par'),'nrScans');
    params.nrScansT2 = readpar_Kea(strcat(parfilestem,'.par'),'nrScansT2');
    params.rxPhase = readpar_Kea(strcat(parfilestem,'.par'),'rxPhase');
    params.rxGain = readpar_Kea(strcat(parfilestem,'.par'),'rxGain');
    params.nrPtsT1 = readpar_Kea(strcat(parfilestem,'.par'),'nrPntsT1');
    params.t1Max = readpar_Kea(strcat(parfilestem,'.par'),'tMax');
    params.t1Est = readpar_Kea(strcat(parfilestem,'.par'),'t1Est');
    params.repTimeT1 = readpar_Kea(strcat(parfilestem,'.par'),'repTime');
    params.repTimeT2 = readpar_Kea(strcat(parfilestem,'.par'),'repTimeT2');
    params.b1Freq = readpar_Kea(strcat(parfilestem,'.par'),'b1Freq');
    params.stepSize = readpar_Kea(strcat(parfilestem,'.par'),'stepSize');
    params.finalDepth = readpar_Kea(strcat(parfilestem,'.par'),'finalDepth');
    params.initDepth = readpar_Kea(strcat(parfilestem,'.par'),'initDepth');
    params.nrEchoesT2 = readpar_Kea(strcat(parfilestem,'.par'),'nrEchoesT2');
    params.nrEchoes = readpar_Kea(strcat(parfilestem,'.par'),'nrEchoes');
    params.smartScan = readpar_Kea(strcat(parfilestem,'.par'),'smartScan');
    
    % Load data
    try
        depthInfo = dlmread('depthInfo.dat');
        depth = depthInfo(:,2);
        T2catch = depthInfo(:,3);
        loadIf = depthInfo(:,4);
        noRough = 0;
    catch EM
        EM
        noRough = 1;
        try
            depth = dlmread('depths.dat');
        catch
            if params.finalDepth == params.initDepth
                depth = params.finalDepth;
            else
                error('proc:BadDepthInfo','The depth info provided by Prospa is incorrect.')
            end
        end
        loadIf = ones(length(depth),1);
    end
    
    %%% Load T1 and T2 data for detailed measurements
    m=1;
    
    if max(loadIf ~= 0)
        
        for n = 1:1:size(loadIf)
            if loadIf(n) == 1
                v = genvarname(['T1data' int2str(depth(n))]);
                w = genvarname(['T2data' int2str(depth(n))]);
                eval([v ' = dlmread(''T1data-' num2str(depth(n)) '.dat'');']);
                eval([w ' = dlmread(''T2data-' num2str(depth(n)) '.dat'');']);
                measDepth(m) = depth(n);
                m = m+1;
            end
        end
        
        allT1Re = zeros(params.nrPtsT1,size(measDepth,2));
        allT1Im = zeros(params.nrPtsT1,size(measDepth,2));
        allT1Mag = zeros(params.nrPtsT1,size(measDepth,2));
        allT2Re = zeros(params.nrEchoesT2,size(measDepth,2));
        allT2Im = zeros(params.nrEchoesT2,size(measDepth,2));
        
        for n = 1:1:size(measDepth,2)
            T1time = eval(['T1data' int2str(measDepth(n)) '(:,1);']);
            T2time = eval(['T2data' int2str(measDepth(n)) '(:,1);']);
            allT1Re(:,n) = eval(['T1data' int2str(measDepth(n)) '(:,2);']);
            allT1Im(:,n) = eval(['T1data' int2str(measDepth(n)) '(:,3);']);
            allT1Mag(:,n) = eval(['sqrt(T1data' int2str(measDepth(n)) '(:,3).^2 + T1data' int2str(measDepth(n)) '(:,2).^2);']);
            allT2Re(:,n) = eval(['T2data' int2str(measDepth(n)) '(:,2);']);
            allT2Im(:,n) = eval(['T2data' int2str(measDepth(n)) '(:,3);']);
        end
        
    end
    
    if autoPhase == 1
        allT2Phased = autophase(complex(allT2Re,allT2Im),1);
        allT1Phased = autophase(complex(allT1Re,allT1Im),1);
    else
        allT1Phased = 0;
        allT2Phased = 0;
    end
    
    if (noRough ~= 1 && params.smartScan == 1);
        for n = 1:1:size(depth,1)
            w = genvarname(['T2data_rough_' int2str(depth(n))]);
            if ispc == 1
                eval([w ' = dlmread(''rough\T2data_rough-' num2str(depth(n)) '.dat'');']);
            else
                eval([w ' = dlmread(''./rough/T2data_rough-' num2str(depth(n)) '.dat'');']);
            end
            T2rough(:,n) = eval(['complex(T2data_rough_' int2str(depth(n)) '(:,2),T2data_rough_' int2str(depth(n)) '(:,3));']);
            if autoPhase == 1
                [T2roughPhased(:,n),T2phaseAngle(:,n)] = autophase(T2rough(:,n),0.1);
            else
                T2roughPhased = 0;
                T2phaseAngle = 0;
            end
        end
    elseif (noRough ~= 1 && params.smartScan ~= 1);
        for n = 1:1:size(depth,1)
            T2rough(:,n) = eval(['complex(T2data' int2str(depth(n)) '(:,2),T2data' int2str(depth(n)) '(:,3));']);
            if autoPhase == 1
                [T2roughPhased(:,n),~] = autophase(T2rough(:,n),0.1);
            else
                T2roughPhased = 0;
                T2phaseAngle = 0;
            end
        end
    end
    
    if m == 1;
        T2time = eval(['T2data_rough_' int2str(depth(1)) '(:,1);']);
    end
    cd(startDir)
end
    