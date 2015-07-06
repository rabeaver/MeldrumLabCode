 function p=satprofile(spectrum,parname)   
    satfreq = readpar(parname,'satfreq'); %number of satcycles
    delaytime = readpar(parname,'ipdelay'); %inter-pulse delay
    pulsetime = readpar(parname,'p1 '); %pulse time; be sure to include space after "p1"
    satcycle = readpar(parname,'satcycle'); %number of satcycles
    temperature = readpar(parname,'temp'); %experimental temperature
    sattime = (pulsetime + delaytime) .* satcycle; %total saturation time including ipdelays
    confidenceInterval = 0.9;

    [C,I] = max(data3);
    data4=squeeze(sum(data3((I(1)-250):(I(1) + 250),:))); %integrate the peak area +/- 250 pts from the highest point
    data4 = data4';
    
    %% Plot of Contrast
    
    for i=2:2:np
        contrast(i/2) = 1- (data4(i-1) - data4(i))/(data4(i-1));
        time(i/2) = sattime(i/2);
        offres(i/2) = data4(i - 1);
        onres(i/2) = data4(i);
    end

    figure
    plot(time,contrast);
    
    
    %% Fitting
    tauon(j,:) = fittingRoutineNoSave(onres,time,confidenceInterval);
    tauoff(j,:) = fittingRoutineNoSave(offres,time,confidenceInterval);
 end