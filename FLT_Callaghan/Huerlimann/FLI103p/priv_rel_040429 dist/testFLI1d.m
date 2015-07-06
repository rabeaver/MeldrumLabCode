

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%				TEST OF FLI1d
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


     clear                                                       % Clear all variables
     

    [data,pp] = readArchmanfile('exp30.ar');
    Tau_1 = imag(data(:,1));

    inputdata = real(data);
	%%%%%%%%%%%%%% Regularization %%%%%%%%%%%%
    Alpha = 0.01;
%	Alpha_Auto = 	0: fixed
%					1: alpha is chosen by BRD method
%					2: a full span of alpha is tested.

    %%%%%%%%%%%%%% Range of T1 and T2 %%%%%%%%%%%%
      
     InitTime_T1 = .001;                      % The initial value of T1 (in seconds)
     FinalTime_T1 = 1;                      % The final value of T1 (in seconds)
     Number_T1 = 100;                        % The number of T1's
            % T1 and T2 are logarithmically spaced
     T1 = logspace(log10(InitTime_T1), log10(FinalTime_T1), Number_T1);

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
        
 		%AllowDCOffset = 1;	% allow a dc offset
        AllowDCOffset = 0; 	% no offset

     	K_1 = Kernel_1 (Tau_1,T1);
     	     	
        [U1, S1, V1] = svds(K_1, 12);
        
  for ijk = 1:size(inputdata,2)
      
            ijk
        Data =  inputdata(:,ijk)';
       
        % use fixed alpha
        [FEst,Alpha_heel,Fit] = FLI1d(Data,K_1,0.005,U1,S1,V1);
        % use t1heel method for regularization
        %[FEst,Alpha_heel,Fit] = FLI1d(Data,K_1,-2,U1,S1,V1);

        alpha_all(ijk) = Alpha_heel;
        error_all(ijk) = std(Fit-Data');
        FEst_all (ijk,:) = FEst(1:Number_T1);
        Data_all(ijk,:) = Data;
        Fit_all  (ijk,:) = Fit';
    
end
figure(1)
    
    figure(1)
    subplot(221)
    plot(Tau_1,Data_all','-',Tau_1,Fit_all')

    subplot(222)
    semilogx(T1,mean(FEst_all,1),T1,std(FEst_all,1),T1,fmodel,'r')
    
    axis tight
    
    subplot(223)
    semilogx(T1,FEst_all')

    subplot(224)
semilogy(alpha_all,'o-')


return


