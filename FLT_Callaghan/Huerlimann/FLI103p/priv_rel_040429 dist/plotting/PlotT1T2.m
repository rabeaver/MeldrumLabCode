% Plot the final estimate of F[][].

     	x = T1; y = T2;
        localpath = pwd;
        
	% Surface plot of f(T1, T2)
	subplot(221)
	surface(x, y, FEst(1:length(y), 1:length(x))),shading interp
	axis square
	h = gca;
	set(h, 'XScale', 'log', 'YScale', 'log')
	xlabel('T1 (secs)', 'FontSize', 9)
	ylabel('T2 (secs)', 'FontSize', 9)
	set(gca, 'FontSize', 9);
	set(gca, 'XTickMode', 'Manual');
	set(gca, 'XTick', [1e-3 1e-2 1e-1 1 10]);
	set(gca, 'YTickMode', 'Manual');
	set(gca, 'YTick', [1e-3 1e-2 1e-1 1 10]);
	title([localpath '/' DataFile], 'FontSize', 9);
	
	v(1) = min(x); v(2) = max(x); V(3) = min(y); v(4) = max(y);
	axis(v)
	h = line([min(x) max(x)], [min(x) max(x)]);
	set(h, 'Color', 'w', 'LineStyle', '-.')
	
    % plot error
    if exist('Fitdata') == 0
        % Compute the best fit
        Fitdata = K_2*FEst*K_1';
    end
    
    subplot(222)
    semilogy(Tau_2,std(Data-Fitdata,0,2))
    hold on
    semilogy(Tau_2,Data(:,end))
    semilogy(Tau_2,mean(abs(Data-Fitdata),2),'r');
    ylabel('std and mean of the fit error along Tau2.')
    title(['Te=' num2str(Tau_2(2)-Tau_2(1)) ',Necho=' num2str(Number_Tau_2)])
    axis tight
	hold off
	
   % Plot T1 and T2 distributions

	subplot(223)
	semilogx(x, T1_dist(1:length(x)), 'or-')
	xlabel('T_1,T_2 (secs)', 'FontSize', 9)
	ylabel(' projections T1(red,o) and T2 (blue) from 2D spectrum,', 'FontSize', 9)
	title(['Tau1=' num2str([Tau_1(1) Tau_1(end)]) 'in ' num2str(Number_Tau_1) ' steps.'], 'FontSize', 9)
	set(gca, 'FontSize', 9);
	hold on
        semilogx(T2,T2_dist, 'b-')
        axis tight
	hold off
	
   % Plot T1/T2 ratio distributions
   	subplot(224)
        if	~exist('T1T2rdist')
            [T1T2rdist,T1T2r] = CalcT1T2r(T1,T2,FEst);
        end
        
        semilogx(T1T2r,T1T2rdist,'o-');
        axis([0.5 100 0 max(T1T2rdist)])
 	xlabel('T1/T2 ratio', 'FontSize', 9)
	ylabel('T1/T2 ratio F(T1/T2) from 2D spectrum', 'FontSize', 9)
	title([date 'T1/T2 distribution with \alpha = ', num2str(Alpha_heel)], 'FontSize', 9)
	set(gca, 'FontSize', 9);

    clear x;
    clear y;
    clear localpath;
    
	orient tall
        
        disp 'to print from matlab, type : print -Pcolor540 -dpsc2';
        
