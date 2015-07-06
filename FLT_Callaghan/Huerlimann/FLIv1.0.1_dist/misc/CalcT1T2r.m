function [dist, ratio] = calcT1T2r(x,y, spec)
% [dist, ratio] = calcT1T2r(x,y, spec)
% T1 along row, T2 along col, spec is NT2 by NT1
% T1 and T2 are logorithmically spaced.

	dlogy = log10(y(2)) - log10(y(1)); dlogx = log10(x(2)) - log10(x(1));
	dd = dlogy*dlogx;
	
	ratio = logspace(-1,2,100);
	dist = zeros(length(ratio),1);
	waste = 0;

	for i = 1:length(y)
		for j = 1:length(x)
			index = (log10(x(j)/y(i)) - log10(ratio(1)) ) * length(ratio) /(log10(ratio(end)/ratio(1)));
			index = round(index);
			if (index < 1) | (index > length(ratio))
				waste = waste +  spec(i,j)*dd;
			else
				
				dist(index) = dist(index) + spec(i,j)*dd;
			end
		
		end
	end

	waste / sum(sum(spec))
