%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                   %
%   Compute marginal distributions, porosity and beta (cos(tipping angle))
%	correction   %
%   for the mis-set 180 pulses                                      %
%                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	
function [x_dist, y_dist, por, beta] = ComputeProjections(FEst,x,y)

length_x = length(x);
length_y = length(y);

dlogx = log10(x(2)/x(1)); 
dlogy = log10(y(2)/y(1));
	
F = FEst(1:length_y, 1:length_x); % T1-T2 distribution

x_dist = sum(F,1) * dlogy; % T1 dist without correction factor

y_dist = sum(F',1) * dlogx; % T2- distribution

if(nargout > 2)
	por = sum(sum(F)) * dlogx * dlogy; % porosity
end

if( nargout > 3)
	beta = sum(FEst(:,end))* dlogx * dlogy/por; % correction factor
end
