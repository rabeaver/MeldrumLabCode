clear
close all
clc

fid=fopen('/Users/tyler/Desktop/test.txt');
fileID = textscan(fid,'%s\t%s\t%u\t%s\t%u\t%s\t%s\t%u\t%u','Delimiter','\t');
fclose(fid);

saveDir = '/Users/tyler/Dropbox/PaintingNewData/PosterFigures/';
%saveFile = 'PaintingList.txt';

%% Process all depthInfo.dat Data and Fit
for masterIndex = 1:length(fileID{1});
% for masterIndex = 1:9;
    clearvars -except masterIndex fileID saveFile saveDir
    close all
    clc
    
    painting_artist = char(fileID{1}(masterIndex));
    painting_title = char(fileID{2}(masterIndex));
    painting_year = fileID{3}(masterIndex);
    painting_color = char(fileID{4}(masterIndex));
    painting_id = fileID{5}(masterIndex);
    parfilestem = char(fileID{6}(masterIndex));
    dataDir = char(fileID{7}(masterIndex));
    depthStart = fileID{8}(masterIndex);
    depthEnd = fileID{9}(masterIndex);
        
    biexp_T2fit = 1;
    T2rough_fit = 1;
    omit_point = 1;
    autoPhase = 1;
    plotRaw = 1;
    saveFigs = 1;
    latexTableOutput = 1;
    
    
    %Load all data
    [params,measDepth,depth,T1time,allT1Re,allT1Im,allT1Mag,T2time,allT2Re,allT2Im,allT1Phased,allT2Phased,T2rough,T2roughPhased] = loadDataFromDepthInfo(dataDir,parfilestem,autoPhase);
     
    
    % Plot all data
    if plotRaw == 1;
%         figure(1)
%         hold all
%         h(1) = subplot(2,2,1);
%         meshc(measDepth,T1time,allT1Re)
%         xlabel('Depth (um)')
%         ylabel('time (ms)')
%         zlabel('signal intensity (arb)')
%         title('T1 Real')
%         h(2) = subplot(2,2,2);
%         meshc(measDepth,T1time,allT1Im)
%         xlabel('Depth (um)')
%         ylabel('time (ms)')
%         zlabel('signal intensity (arb)')
%         title('T1 Imag')
%         h(3) = subplot(2,2,3);
%         meshc(measDepth,T2time(omit_point+1:end),allT2Re(omit_point+1:end,:))
%         xlabel('Depth (um)')
%         ylabel('time (ms)')
%         zlabel('signal intensity (arb)')
%         title('T2 Real')
%         h(4) = subplot(2,2,4);
%         meshc(measDepth,T1time,allT1Mag)
%         xlabel('Depth (um)')
%         ylabel('time (ms)')
%         zlabel('signal intensity (arb)')
%         title('T1 Mag')
% 
%       

    % Detailed T2 plot
    profileFig = figure('Units','normalized','Position',[0.66 0.5 0.33 0.45]);   

    hold on
    plot(depth,real(T2roughPhased(omit_point+1,:)),'-k')
    line([min(depth);max(depth)],[0;0],'Color','k')
%     ylim([-1e4 8e4]);
    
    for i=1:length(measDepth)
        xind = measDepth(i);
        yind = real(T2roughPhased(omit_point+1,find(depth == measDepth(i))));
        indlabel = num2str(i);
        
        try
            yprev = real(T2roughPhased(omit_point+1,find(depth == measDepth(i))-1));
        catch
            yprev = yind;
        end
        
        try
            yafter = real(T2roughPhased(omit_point+1,find(depth == measDepth(i))+1));
        catch
            yafter = yind;
        end
        
        
        if(yind <= yprev && yind <= yafter) %min
            text(xind,yind,indlabel,'VerticalAlignment','top','HorizontalAlignment','center');
        elseif(yind <= yprev && yind > yafter) %pos slope
            text(xind,yind,indlabel,'VerticalAlignment','middle','HorizontalAlignment','right');
        elseif(yind > yprev && yind <= yafter) %neg slope
            text(xind,yind,indlabel,'VerticalAlignment','middle','HorizontalAlignment','left');
        elseif(yind > yprev && yind > yafter) %max
            text(xind,yind,indlabel,'VerticalAlignment','bottom','HorizontalAlignment','center');
        end
    end
        

%         figure(3)
%         mesh(depth,T2time(omit_point+1:end),real(T2roughPhased(omit_point+1:end,:)))
    end
    
    
    % Fit rough data to T2 decay    
    if T2rough_fit == 1       
        guesses = [0;max(max(real(T2roughPhased)));T2time(round(length(T2time)/10))];
            CI = 90; %desired confidence interval in percent
            ypred_T2_rough = zeros(1001,length(depth));
            y0_T2_rough = zeros(length(depth),2);
            A_T2_rough = zeros(length(depth),2);
            tau_T2_rough = zeros(length(depth),2);
            
            for n = depthStart:depthEnd; %1:length(depth)
                try
                    [xfit_T2,ypred_T2_rough(:,n),T2_coeffs] = monodecay_t2fit(T2time(omit_point+1:end),real(T2roughPhased(omit_point+1:end,n)),guesses,CI);
                    y0_T2_rough(n,:) = T2_coeffs(1,1:2);
                    A_T2_rough(n,:) = T2_coeffs(2,1:2);
                    tau_T2_rough(n,:) = T2_coeffs(3,1:2);
                catch 
                    y0_T2_rough(n,:) = [0;0];
                    A_T2_rough(n,:) = [0;0];
                    tau_T2_rough(n,:) = [0;0];
                end
            end          
        tau_T2_rough(:,3) = tau_T2_rough(:,2)./tau_T2_rough(:,1)*100;
    end
   
    tau0profile = y0_T2_rough(:,1) + A_T2_rough(:,1);
    tau30profile = y0_T2_rough(:,1) + A_T2_rough(:,1).*exp(-0.03./tau_T2_rough(:,1));
    tau100profile = y0_T2_rough(:,1) + A_T2_rough(:,1).*exp(-0.1./tau_T2_rough(:,1));
    
    profileExtrap = figure('Units','normalized','Position',[0.66 0 0.33 0.45]); 
    hold on
    plot(depth,tau0profile)
    plot(depth,tau30profile,'-k')
    plot(depth,tau100profile,'-r')
    
        for i=1:length(measDepth)
        xind = measDepth(i);
        yind = tau0profile(find(depth == measDepth(i)));
        indlabel = num2str(i);
        
        try
            yprev = real(T2roughPhased(omit_point+1,find(depth == measDepth(i))-1));
        catch
            yprev = yind;
        end
        
        try
            yafter = real(T2roughPhased(omit_point+1,find(depth == measDepth(i))+1));
        catch
            yafter = yind;
        end
        
        
        if(yind <= yprev && yind <= yafter) %min
            text(xind,yind,indlabel,'VerticalAlignment','top','HorizontalAlignment','center');
        elseif(yind <= yprev && yind > yafter) %pos slope
            text(xind,yind,indlabel,'VerticalAlignment','middle','HorizontalAlignment','right');
        elseif(yind > yprev && yind <= yafter) %neg slope
            text(xind,yind,indlabel,'VerticalAlignment','middle','HorizontalAlignment','left');
        elseif(yind > yprev && yind > yafter) %max
            text(xind,yind,indlabel,'VerticalAlignment','bottom','HorizontalAlignment','center');
        end
        end
    
        
    
%     break
    
% %     testdata = y0_T2_rough(:,1) + A_T2_rough(:,1);
% %     testratio = testdata./real(T2roughPhased(2,:))';
% %     figure(profileFig)
% %     hold on
% %     plot(depth,testdata,'-ob')
% % %     ylim([-1e3 10e4]);
% %     
% %     figure
% %     plot(depth,testratio,'-or')
% % %     ylim([0 20])
    
    
    % Fit regular data to T2 decay   
    guesses = [0;max(max(allT2Re));T2time(round(length(T2time)/8))];  
    CI = 90; %desired confidence interval in percent  
    ypred_T2 = zeros(1001,size(measDepth,2));
    
    for n = 1:1:size(measDepth,2)
        try
            [xfit_T2,ypred_T2(:,n),T2_coeffs,T2_resid(:,n)] = monodecay_t2fit(T2time(omit_point+1:end),real(allT2Phased(omit_point+1:end,n)),guesses,CI);
            y0_T2(n,:) = T2_coeffs(1,1:2);
            A_T2(n,:) = T2_coeffs(2,1:2);
            tau_T2(n,:) = T2_coeffs(3,1:2);
        catch
            y0_T2(n,:) = [0,0];
            A_T2(n,:) = [0,0];
            tau_T2(n,:) = [0,0];
            T2_resid(:,n) = zeros(length(T2time)-omit_point,1);
        end
    end
    
    tau_T2(:,3) = tau_T2(:,2)./tau_T2(:,1)*100; 
    
    % optional biexponential fit    
    if biexp_T2fit == 1; 
        guesses = [0;max(max(allT2Re));T2time(round(length(T2time)/20));max(max(allT2Re))/4;T2time(round(length(T2time)/5))];
        CI = 90; %desired confidence interval in percent
        ypred_T2_bi = zeros(1001,size(measDepth,2));
        
        for n = 1:1:size(measDepth,2)
            try
                [xfit_T2,ypred_T2_bi(:,n),T2_coeffs_bi,T2_resid_bi(:,n)] = bidecay_t2fit(T2time(omit_point+1:end),real(allT2Phased(omit_point+1:end,n)),guesses,CI);
                y0_T2_bi(n,:) = T2_coeffs_bi(1,1:2);
                A_T2_bi(n,:) = T2_coeffs_bi(2,1:2);
                tau_T2_bi(n,:) = T2_coeffs_bi(3,1:2);
                A2_T2_bi(n,:) = T2_coeffs_bi(4,1:2);
                tau2_T2_bi(n,:) = T2_coeffs_bi(5,1:2);
            catch
                guesses = [y0_T2(n,1);A_T2(n,1);tau_T2(n,1);A_T2(n,1)/5;tau_T2(n,1)*5];
                try
                    [xfit_T2,ypred_T2_bi(:,n),T2_coeffs_bi,T2_resid_bi(:,n)] = bidecay_t2fit(T2time(omit_point+1:end),real(allT2Phased(omit_point+1:end,n)),guesses,CI);
                    y0_T2_bi(n,:) = T2_coeffs_bi(1,1:2);
                    A_T2_bi(n,:) = T2_coeffs_bi(2,1:2);
                    tau_T2_bi(n,:) = T2_coeffs_bi(3,1:2);
                    A2_T2_bi(n,:) = T2_coeffs_bi(4,1:2);
                    tau2_T2_bi(n,:) = T2_coeffs_bi(5,1:2);
                catch
                    y0_T2_bi(n,:) = [0,0];
                    A_T2_bi(n,:) = [0,0];
                    tau_T2_bi(n,:) = [0,0];
                    A2_T2_bi(n,:) = [0,0];
                    tau2_T2_bi(n,:) = [0,0];
                end
            end
            
        end
            
            tau_T2_bi(:,3) = tau_T2_bi(:,2)./tau_T2_bi(:,1)*100;
            tau2_T2_bi(:,3) = tau2_T2_bi(:,2)./tau2_T2_bi(:,1)*100;
      
    end
    
   
    
    % Fit to T1 decay   
    guesses = [0;max(max(real(allT1Phased)));T1time(round(length(T1time)/6))];
    CI = 90; %desired confidence interval in percent
    ypred_T1 = zeros(1001,size(measDepth,2));
    
    for n = 1:1:size(measDepth,2)
        [xfit_T1,ypred_T1(:,n),T1_coeffs,T1_resid(:,n)] = T1_fit(T1time,real(allT1Phased(:,n)),guesses,CI);
        y0_T1(n,:) = T1_coeffs(1,1:2);
        A_T1(n,:) = T1_coeffs(2,1:2);
        tau_T1(n,:) = T1_coeffs(3,1:2);
    end
    
    tau_T1(:,3) = tau_T1(:,2)./tau_T1(:,1)*100;
    
    % Plot data and fits and select which to save
    
    T1T2fig = figure('Units','normalized','Position',[0 0 0.66 0.66]); %'Units','normalized','OuterPosition',[0.15 0.15 0.7 0.7]); %,'MenuBar','none');
    n = size(measDepth,2);
    str{1} = '0';
    
    for i = 1:n
        str{i} = num2str(i);
        subplot(6,n,[i i+n])
        hold on
        %     plot(xfit_T1,ypred_T1(:,i))
        %     plot(T1time,allT1Re(:,i),'ob')
        %     plot(T1time,allT1Im(:,i),'or')
        %     plot(T1time,allT1Mag(:,i),'ob')
        plot(T1time,real(allT1Phased(:,i)),'ok')
        plot(xfit_T1,ypred_T1(:,i),'-k')
        ylim([-0.5*max(max(allT1Mag)) 1.1*max(max(allT1Mag))])
        xlim([0 max(T1time)])
        title(i)
        
        subplot(6,n,2*n+i)
        hold on
        plot(T1time,T1_resid(:,i),'-ok')
        line([0;max(T1time)],[0;0],'Color','k')
        ylim([-1.1*max(max(T1_resid)) 1.1*max(max(T1_resid))])
        xlim([0 max(T1time)])
        
        
        subplot(6,n,[3*n+i i+4*n])
        hold on
        %     plot(T2time,allT2Re(:,i),'-ob')
        %     plot(T2time,allT2Im(:,i),'-or')
        plot(T2time,real(allT2Phased(:,i)),'ok')
        plot(xfit_T2,ypred_T2(:,i),'-k')
        if biexp_T2fit == 1;
            plot(xfit_T2,ypred_T2_bi(:,i),'-.k')
        end
        %     plot(xfit_T2_lin,ypred_T2_lin(:,i),'-r')
        limits = max(max(real(allT2Phased)));
        ylim([-0.1*limits 1.1*limits])
        xlim([0 max(T2time)])
        
        subplot(6,n,5*n+i)
        hold on
        plot(T2time(omit_point + 1:end),T2_resid(:,i),'-k')
        plot(T2time(omit_point + 1:end),T2_resid_bi(:,i),'-.k')
        %     plot(T2time(omit_point + abs(linearizedFit-1):end),T2_resid(:,i),'-b')
        line([0;max(T2time)],[0;0],'Color','k')
        %     bar([T2time(omit_point+1:end),T2_resid(:,i),T2_resid_bi(:,i)])
        %     bar(T2time(omit_point+1:end),T2_resid_bi(:,i))
        ylim([-0.2*limits 0.2*limits])
        xlim([0 max(T2time)])
    end
   
%     break
    
    % List dialog box to select which datasets to use
    [s,saveData] = listdlg('PromptString','Save datasets:','OKString','Save',...
        'ListString',str,'ListSize',[140 100]);
    
    if saveFigs == 1 && saveData == 1;
        filename_prof = strcat(saveDir,'Profile_',painting_artist,'_',painting_color,'_',num2str(painting_id),'.eps');
        filename_profExtrap = strcat(saveDir,'Profile_Extrapolated_',painting_artist,'_',painting_color,'_',num2str(painting_id),'.eps');
        filename_t1t2 = strcat(saveDir,'T1T2_',painting_artist,'_',painting_color,'_',num2str(painting_id),'.eps');
        export_fig(filename_prof,profileFig)
        export_fig(filename_profExtrap,profileExtrap)
        export_fig(filename_t1t2,T1T2fig)
    end
        
    close all
    
    % Save data to file
    filename = strcat(saveDir,'text_',painting_artist,'_',painting_color,'_',num2str(painting_id),'.txt');
    latexFilename = strcat(saveDir,'table_',painting_artist,'_',painting_color,'_',num2str(painting_id),'.tex');
    if biexp_T2fit ~= 1;
        outputCoeffs = [y0_T1, A_T1, tau_T1, y0_T2, A_T2, tau_T2];
        if saveData == 1;
            fid = fopen(filename, 'a');
            for i = 1:length(s)
                fprintf(fid, '%s\t %s\t %4i\t %s\t %u\t %u\t %6f\t %6f\t %6f\t %6f\t %6f\t %6f\t %6f\t %6f\t %6f\t %6f\t %6f\t %6f\t %6f\t %6f\t\n', painting_title, painting_artist, painting_year, painting_color, painting_id, s(i), outputCoeffs(s(i),:));
            end
            fclose(fid);
        end
    else
        outputCoeffs = [y0_T1, A_T1, tau_T1, y0_T2, A_T2, tau_T2,y0_T2_bi, A_T2_bi, tau_T2_bi, A2_T2_bi, tau2_T2_bi];
        if saveData == 1;
            fid = fopen(filename, 'a');
            for i = 1:length(s)
                fprintf(fid, '%s\t %s\t %4i\t %s\t %u\t %u\t %6f\t %6f\t %6f\t %6f\t %6f\t %6f\t %6f\t %6f\t %6f\t %6f\t %6f\t %6f\t %6f\t %6f\t %6f\t %6f\t %6f\t %6f\t %6f\t %6f\t %6f\t %6f\t %6f\t %6f\t %6f\t %6f\t\n', painting_title, painting_artist, painting_year, painting_color, painting_id, s(i), outputCoeffs(s(i),:));
            end
            fclose(fid);
        end
    end
    
    if latexTableOutput == 1 && saveData == 1;
        precis = 5; %precision for the num2str output
        newstr = 'cc';
        str = 'c';
        headerstr = ' & ';
        t1_y0str = ' ';
        t1_Astr = ' ';
        t1_taustr = ' ';
        t2_y0str = ' ';
        t2_Astr = ' ';
        t2_taustr = ' ';
        t2bi_y0str = ' ';
        t2bi_A1str = ' ';
        t2bi_tau1str = ' ';        
        t2bi_A2str = ' ';
        t2bi_tau2str = ' ';        
        
        for i=1:length(s);
            headeradd = [' & ',num2str(s(i)),' '];
            newstr = strcat(newstr,str);
            headerstr = strcat(headerstr,headeradd);
            t1_y0str = strcat(t1_y0str,' & ',num2str(outputCoeffs(s(i),1),precis),' ');
            t1_Astr = strcat(t1_Astr,' & ',num2str(outputCoeffs(s(i),3),precis),' ');
            t1_taustr = strcat(t1_taustr,' & ',num2str(outputCoeffs(s(i),5),precis),' ');            
            t2_y0str = strcat(t2_y0str,' & ',num2str(outputCoeffs(s(i),8),precis),' ');
            t2_Astr = strcat(t2_Astr,' & ',num2str(outputCoeffs(s(i),10),precis),' ');
            t2_taustr = strcat(t2_taustr,' & ',num2str(outputCoeffs(s(i),12),precis),' ');  
            t2bi_y0str = strcat(t2bi_y0str,' & ',num2str(outputCoeffs(s(i),15),precis),' ');
            t2bi_A1str = strcat(t2bi_A1str,' & ',num2str(outputCoeffs(s(i),17),precis),' ');
            t2bi_tau1str = strcat(t2bi_tau1str,' & ',num2str(outputCoeffs(s(i),19),precis),' ');     
            t2bi_A2str = strcat(t2bi_A2str,' & ',num2str(outputCoeffs(s(i),22),precis),' ');
            t2bi_tau2str = strcat(t2bi_tau2str,' & ',num2str(outputCoeffs(s(i),24),precis),' ');                             
        end
        
        latexString = ['\\begin{tabular}{%s} \\\\ \n',... %newstr, sets up columns
                       '\\toprule\n',...
                       ' & & \\multicolumn{%u}{c}{Peak number} \\\\ \n',... %number of columns to center title over
                       '%s \\\\ \n',... %headerstr, sets up peak numbers with pm signs
                       '\\midrule\n',...                      
                       '\\multirow{3}{*}{$T_1$}',' & $y_0$ ',t1_y0str,'\\\\ \n',...
                       ' & $A$ ',t1_Astr,'\\\\ \n'...
                       ' & $\\tau$ ',t1_taustr,'\\\\ \n'...                       
                       '\\midrule\n',...                      
                       '\\multirow{3}{*}{$T_2$}',' & $y_0$ ',t2_y0str,'\\\\ \n',...
                       ' & $A$ ',t2_Astr,'\\\\ \n'...
                       ' & $\\tau$ ',t2_taustr,'\\\\ \n'...  
                       '\\midrule\n',...                      
                       '\\multirow{5}{*}{$T_2$ bi}',' & $y_0$ ',t2bi_y0str,'\\\\ \n',...
                       ' & $A_1$ ',t2bi_A1str,'\\\\ \n'...
                       ' & $\\tau_1$ ',t2bi_tau1str,'\\\\ \n'... 
                       ' & $A_2$ ',t2bi_A2str,'\\\\ \n'...
                       ' & $\\tau_2$ ',t2bi_tau2str,'\\\\ \n'...                        
                       '\\bottomrule\n',...
                       '\\end{tabular}'];
                   
        fid = fopen(latexFilename, 'w');
            fprintf(fid, latexString, newstr, length(s), headerstr);
        fclose(fid);
    end
end
