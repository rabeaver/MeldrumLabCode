close all
clc
clear



for ll = 1:5
    filedir = '/Users/tyler/Desktop/Hamada Samples/Profiles/Camphene/Second/';
    filename = 'Camphene with Rubber.2d';
    parloc = strcat(filedir,num2str(ll),'/acqu.par');
    
    par.echoTime = readpar_Kea(parloc,'echoTime');
    par.stepSize = readpar_Kea(parloc,'stepSize');
    par.initDepth = readpar_Kea(parloc,'initDepth');
    par.finalDepth = readpar_Kea(parloc,'finalDepth');
    
    depth = par.initDepth:par.stepSize:par.finalDepth; %um
    
    
    data = readKea4d(strcat(filedir,num2str(ll),'/',filename));
    dd = data.data;
    dd = reshape(dd,data.xDim,data.yDim);
    echoVec(:,ll) = 1e-3*(par.echoTime:par.echoTime:data.xDim*par.echoTime); %ms
    
    for ii = 1:data.xDim
        ee(ii,:,ll) = smooth(dd(ii,:),7);
    end

    
    figure(1)
    hold on
    plot(depth,ee(2,:,ll))
    set(gca,'XDir','Reverse')
    xlim([par.finalDepth,par.initDepth]);
    legend('1','2','3','4','5','6')
    
    
    
end



