% Plot the final estimate of F[][].


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                    %
%   Compute and plot T1 and T2 distributions  and beta (correction for mis-set 180   %
%                                                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  PlotFEstimate(x, y, FEst, alpha,DataName)  

     	
	% Surface plot of f(T1, T2)
	subplot(211)
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
	title([DataName, ':F(T1, T2) for \alpha  = ', num2str(alpha)], 'FontSize', 9);
	
	v(1) = min(x); v(2) = max(x); V(3) = min(y); v(4) = max(y);
	axis(v)
	h = line([min(x) max(x)], [min(x) max(x)]);
	set(h, 'Color', 'w', 'LineStyle', '-.')
	
%	% Compute T1, T2 distributions, porosity and beta
	
	[x_dist, y_dist, por, beta] = ComputeProjections(FEst, x,y);
	
   % Plot T1 and T2 distributions
	%figure('Units', 'Inches','Position', [0.5 0.5 7.5 10])
	subplot(223)
	semilogx(x, x_dist(1:length(x)), 'r-.')
	xlabel('T_1 (secs)', 'FontSize', 9)
	ylabel(' F(T_1)', 'FontSize', 9)
	title(['T_1 distribution with \alpha = ', num2str(alpha)], 'FontSize', 9)
	set(gca, 'FontSize', 9);
	
	subplot(224)
	semilogx(y, y_dist, 'r-.')
	xlabel('T_2 (secs)', 'FontSize', 9)
	ylabel('F(T_2)', 'FontSize', 9)
	title(['T_2 distribution with \alpha = ', num2str(alpha)], 'FontSize', 9)
	set(gca, 'FontSize', 9);
	fprintf(1, ' por. = %g\n',por);

	orient tall
