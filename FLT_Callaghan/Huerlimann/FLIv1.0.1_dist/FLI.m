

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                    %
%                                                                                    %
%        Main matlab file to do T1-T2 inversion.  It solves the problem of the       %
%        form M(tau_2, tau_1) = double integral ( F(T1, T2) (1.0 -                   %
%        2exp(-tau_1/T1)) exp(-tau_2/T2)) dT1 dT2                                    %
%                                                                                    %
%        where tau_1 corresponds to the different wait times and tau_2               %
%        corresponds to the equally spaced times at which CPMG data is acquired.     %
%                                                                                    %
%        Here, M(tau_2, tau_1) is measured experimentally by following a CPMG        %
%        following inversion recovery.  The objective is to estimate F(T1, T1)       %
%        subject to the constraint that F be greater than or equal to zero, given    %
%        the data.                                                                   %
%                                                                                    %
%        The problem is a least-squares problem and can be written to be of the      %
%        form, D = K2 * F * K1' + E where D(i,j) has the data at the i-th tau_2      %
%        value and j-th tau_1 value K2(i,j) = exp(-tau_2(i)/T2(j)) F(i,j) is the     %
%        (i,j)th element of F K1(i,j) = 1.0 - 2*exp(-tau_1(j)/T1(i)) E refers to     %
%        additive white Gaussian noise.                                              %
%                                                                                    %
%        The solution is in two steps :                                              %
%                                                                                    %
%        1.  Data compression : The data compression is done using singular value    %
%        decomposition of K2 and K1.  Consider r significant singular values of      %
%        K2 and p signifianct singular values of K1.  It can be shown that (a)       %
%        the data can be compressed to be of the size (r*p).  Let the compressed     %
%        data be denoted by D_.  (b) the problem can be re-formulated to be          %
%                                                                                    %
%        D_ = K2_ * F * K1'_ - Equation (a)                                          %
%                                                                                    %
%        where K2_ and K1_ are matrices and are found from the SVD of K2 and K1.     %
%                                                                                    %
%        2.  Solve (a) with zero-th order regularization with procedure given in     %
%        Butler, Reeds, and Dawson paper, subject to F>=0.                           %
%                                                                                    %
%        Inputs to the file are :                                                    %    
%        1. Number_T1, Number_T2, InitTime_T1, InitTime_T2, Final_T1, Final_T2,      %
%        2. The number of alphas, alpha_min, alpha_max                               %
%        3. Flag that indicates of data is being simulated or data already exists    %
%           as an experimental data set                                              %
%        4. Flag that indicates if the singular values are being computed or have    %
%           have been pre-computed and saved ; if flag is on, the upper bounds on    %
%           the number of singularvalues required for S1 and S2.                     %
%        5. If simulator is being used : InitTime_Tau_1, Final_Time_Tau_1,           %
%           Number_Tau_1, Echospacing, InitTime_Tau_2, Final_Time_Tau_2,             %
%           Number_Tau_2, Name of file where the data is stored                      %
%        6. The name of the data file (should have data, tau_1, tau_2, noisestd)     %
%        7. Name of file where the SVD is stored (if it is stored)                   % 
%        8. Ensure that if T1 or T2 or tau_1 or tau_2 changed, we have to compute    %
%           the SVD again.                                                           %
%                                                                                    %
%                                                                                    %
%                                                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Original March 2001 
%	modified April 2001, September 2002.
%	
%	Rename for FLI.m for fast Laplace Inversion
%
%	Details were published in 
%	L. Venkataramanan et. al., IEEE Tran. Signal Proc. 50, 1017-1026 (May, 2002).	 %
%	Y.-Q. Song, et al., J. Magn. Reson. 154, 261-268(2002).			    			 %
%	M D Hurlimann and L Venkataramanan, J. Magn. Reson. 157, 31 (2002).		    	 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%				BEGINNING OF THE MAIN PROGRAM
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


     clear                                                       % Clear all variables
     localpath1 = pwd;
     addpath([localpath1 filesep 'plotting']);
     addpath([localpath1 filesep 'misc']);
     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  USER INPUT FILE  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

     % Name of data file : data has to be in a *.mat file with four fields :
     % Tau_1, Tau_2, NoiseStd, Data
     % Result file will be store in the same location.


	R = input('Please provide the input data file name:','s');
	
	if isempty(R)
		DataFile=['examples' filesep 'example1.mat']
	else
		DataFile = R
	end
	clear('R');

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  USER MODIFICATIONS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




	%%%%%%%%%%%%%% Regularization %%%%%%%%%%%%
    Alpha_Auto = 0;
%	Alpha_Auto = 	0: fixed
%					1: alpha is chosen by BRD method
%					2: a full span of alpha is tested.
%
    AlphaStart = 1;                        % Minimum value of log10(Alpha)
%if fixed alpha
%    AlphaStart = 1;                        % Minimum value of log10(Alpha)
%    Alpha_Auto = 0;                % Flag that indicates if alpha is automatically chosen
                                    % 0 : if alpha is fixed
                                    % 1 : if Alpha is selected by the BRD method

      
      
	%%%%%%%%%%%%%% Range of T1 and T2 %%%%%%%%%%%%
      
     InitTime_T1 = .0001;                      % The initial value of T1 (in seconds)
     FinalTime_T1 = 100;                      % The final value of T1 (in seconds)
     Number_T1 = 32;                        % The number of T1's
       
     InitTime_T2 = .0001;                      % The initial value of T2 (in seconds)
     FinalTime_T2 =100;                      % The final value of T2(in seconds)
     Number_T2 = 32;                        % The number of T2's




	%%%%%%%%%%%%% Definitions of Kernel functions %%%%%%%%%%%%%%%
	
	%%%% CONVENTION
		% Tau is a vertical vector, TimeConst is horizontal vector,
	% and the resulting kernels are matrices.
	
	%%%%%%%%%%%% Kernel along the first dimension %%%%%%%%%%%%%%
	%exponential decay
	Kernel_1 = inline('exp(- Tau * (1./ TimeConst))','Tau','TimeConst');
	
	%inv recovery
	%Kernel_1 = inline('1- 2*exp( - Tau * (1./ TimeConst))','Tau','TimeConst');
	
		
		%%%%%%%%%%%% Kernel modification for DC offset %%%%%%%%%%%%%

        % for inversion-recovery data, include the a DC offset in the kernel 
        % will help estimate the error in the pi pulse. However, if pi error is too
        % large, it is better to modify the kernel function to account for it first.
        
        % the kernel will be modified by adding a column of 1 as the last column
        % so that the program will find the pi pulse error, or other effects to produce
        % a dc  offset.
        
 		AllowDCOffset = 1;	% allow a dc offset
       % AllowDCOffset = 0; 	% no offset


        
	%%%%%%%%%%%% Kernel along the second dimension %%%%%%%%%%%%%%
    %exponential decay
	Kernel_2 = inline('exp( - Tau * (1./ TimeConst))','Tau','TimeConst');



       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% END OF USER INPUTS & MODIFICATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








     %  ComputeSingularValues is a flag that indicates if SVD needs to be computed
     %  if flag is turned off, it is assumed that K2, K1 and their SVDs are stored in 
     %  file SingularValues.mat
	 %     ComputeSingularValues = 1;  
     Number_S1 = 10;                % Upper limit for number of singular values for K1
     Number_S2 = 10;                % Upper limit for number of singular values for K2
	
     
     S = load  (DataFile);
     Tau_1 = S.Tau_1; Tau_2 = S.Tau_2; Data = S.Data; NoiseStd = S.NoiseStd;
     clear S
     
     Number_Tau_1 = length(Tau_1);
     Number_Tau_2 = length(Tau_2);
     ConditionNumber = 1000;




     %      Analysis of the data :
     
     fprintf(1,' Beginning analysis of the data ...\n');
     
     % T1 and T2 are logarithmically spaced
     T1 = logspace(log10(InitTime_T1), log10(FinalTime_T1), Number_T1);
     T2 = logspace(log10(InitTime_T2), log10(FinalTime_T2), Number_T2);
          
     ComputeSingularValues = 1;	%% always
     
     if (ComputeSingularValues)    	
     
     	% Compute matrix K_2 where K_2(Tau_2, T2) = exp(-Tau_2/T2) dlogT2  
     	fprintf(1,'Analysis ... Computing K_2 ...\n');

     	K_2 = Kernel_2 (Tau_2,T2);
     	
     	 	
     	fprintf(1,'Analysis ... Computing K_1 ...\n');     	
     	
     	%K_1 = 1 - exp(-Tau_1 * (1 ./T1) ) ;	%inv recovery
     	%K_1 = exp(-Tau_1 * (1 ./T1) ) ;		%exponential decay
     	K_1 = Kernel_1 (Tau_1,T1);
     	
     	
     	% Adjust for imperfect 180 pulse, in an inversion recovery experiment.
        if AllowDCOffset == 1
            for i = 1:Number_Tau_1  K_1(i,Number_T1+1) = 1.0; end
            fprintf(1,'Allow DC offset.. \n');
        else
            for i = 1:Number_Tau_1  K_1(i,Number_T1+1) = 0.0; end
	end
        
     	% SVD of K_1
     	fprintf(1,'Computing SVD of K_1 ... ');
     	t = cputime;
     	[U1, S1, V1] = svds(K_1, Number_S1); 
     	e = cputime -t;
     	fprintf(1, 'Time for svd for K_1 = 	%d\n',e);
     	
     	%SVD of K2
     	fprintf(1,'Computing SVD of K_2 ...');
     	t = cputime;
     	[U2, S2, V2] = svds(K_2, Number_S2); 
     	e = cputime -t;
     	fprintf(1, 'Time for svd for K_2 = %d\n',e);
     	
     	% Plot the singular values
     	%PlotSingularValues(S1, S2);
 
	end  


    switch Alpha_Auto 
        
        case {0,1}	% 0:for a fixed alpha, 
                        % 1:BRD method
             % Compute the best estimate of F, subject to F>=0.
            [FEst, CompressedData, Chi, Alpha_heel] = ...
                FLIEstimate(Data, Tau_1, Tau_2,  U1, U2, V1, V2, ...
                    S1, S2, T1, T2, AlphaStart, NoiseStd, Alpha_Auto, ConditionNumber);
 
            WriteFileName = strcat(DataFile,'Re.mat');			%% result is stored.
            save(WriteFileName);
            fprintf(1,'The matlab environment has been saved in %s.\n',WriteFileName);


           % 	compute projects, T1 dist and T2 dist, 
            %	porosity and tipping angle correction
            [T1_dist, T2_dist, por, beta] = ComputeProjections(FEst, T1, T2);
			fprintf('beta (correction for mis-set 180 pulse) = %d\n', beta);
			
            % Compute the best fit
            % Fitdata has the same dimensions as Data.
            
            Fitdata = K_2*FEst*K_1';
                               
            %if everything is ok, save again
            WriteFileName = strcat(DataFile,'Re.mat');
            save(WriteFileName);
	        fprintf(1,'The matlab environment has been saved in %s.\n',WriteFileName);            

        case 2	 	% full span of alpha
                AlphaSpan = logspace(-3,10,100);	%% calc for a series of Alpha. the range may be changed 
                					%% to suit the specific needs
                AlphaSpan = fliplr(AlphaSpan);
                Alpha_AutoTmp = 0;
                for ijk = 1:length(AlphaSpan)
                    [FEst, CompressedData, Chi, Alpha_heel] = FLIEstimate(Data, Tau_1, Tau_2,  U1, U2, V1, V2, ...
                    S1, S2, T1, T2, AlphaSpan(ijk), NoiseStd, Alpha_AutoTmp, ConditionNumber);
                    CompressedError(ijk) = Chi;
                    FDensity{ijk} = FEst;
                end 
                
                
                % findheel	- FindHeel.m
                CompressedError = CompressedError(1:end-5);
                AlphaSpan = AlphaSpan(1:end-5);
                logChi = fliplr(log10(CompressedError));
                AlphaSpan = fliplr(AlphaSpan);
                logalpha = log10(AlphaSpan(2)/AlphaSpan(1));
                
                for i = length(logChi) : -1 :2
                    dlogChi(i) = (logChi(i) - logChi(i-1))/logalpha;
                end
                dlogChi(1) = 0;
                
                 
                j = find(dlogChi >= 0.1);
                
               Alpha_heel = AlphaSpan(j(1));
                
                [FEst, CompressedData, Chi, Alpha_heel] = ...
                    FLIEstimate(Data, Tau_1, Tau_2,  U1, U2, V1, V2, ...
                        S1, S2, T1, T2,Alpha_heel, NoiseStd, Alpha_AutoTmp, ConditionNumber);
                
                WriteFileName = strcat(DataFile,'ReFull.mat');				%% results stored.
                save(WriteFileName);
                fprintf(1,'The matlab environment has been saved in %s.\n',WriteFileName);

				%%%%%  plot chi-sq
                figure
                plot(log10(AlphaSpan), logChi)
                hold on
                plot(log10(AlphaSpan), logChi, '*')
                
               
                
                 figure
                plot(log10(AlphaSpan), dlogChi)

			   % 	compute projects, T1 dist and T2 dist, 
				%	porosity and tipping angle correction
				[T1_dist, T2_dist, por, beta] = ComputeProjections(FEst, T1, T2);
				fprintf('beta (correction for mis-set 180 pulse) = %d\n', beta);
				
				% Compute the best fit
				% Fitdata has the same dimensions as Data.
				
				Fitdata = K_2*FEst*K_1';

 
                %if everything is ok, save again
                WriteFileName = strcat(DataFile,'ReFull.mat');
                save(WriteFileName);

             otherwise % SWITCH
            	fprintf(1,'Alpha_Auto = %d. Not a valid choice.\n',Alpha_Auto);
        end	%SWITCH
        

		%%%%%%%%%%% PLOTTING %%%%%%%%%%%%%%%%%
        % calculate T1T2 ratio
        [T1T2rdist,T1T2r] = CalcT1T2r(T1,T2,FEst);
        
        
	% Plot f(T1,T2) and projections/marginal distributions.

	% PlotFEstimate(T1, T2, FEst, Alpha_heel,DataFile);
	PlotT1T2
				
	% Plot six sections of the 2-D data to ensure that the fit looks reasonable
	% PlotData2(Tau_1, Tau_2, Data, Fitdata);


	rmpath([localpath1 filesep 'plotting']);
        rmpath([localpath1 filesep 'misc']);
