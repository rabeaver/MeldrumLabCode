
  
function PlotFig(Tau_1, Tau_2, Data, label, titlestring, zmin, zmax, ...
	color, xlog, ylog)

	
if ( (length(Tau_1)  == 1)  & (length(Tau_2) > 1) )
	plot(Tau_2, Data, color)
	set(gca, 'XScale', xlog, 'YScale', ylog)
	xlabel('\tau_2')
	ylabel(label)
	v = axis; v(3) = zmin; v(4) = zmax;
	axis(v)
	title(titlestring)
elseif ( (length(Tau_2) == 1)  & (length(Tau_1) > 1) )
	plot(Tau_1, Data, color)
	set(gca, 'XScale', xlog, 'YScale', ylog)
	xlabel('\tau_1')
	ylabel(label)
	v = axis; v(3) = zmin; v(4) = zmax;
	axis(v)
	title(titlestring)
elseif ( (length(Tau_1) == 1) & (length(Tau_2) == 1) )
	[m,n] = size(Data);
	if  ( (m==1) | (n==1) )
		plot(Data, color)	
	else 
		mesh(Data) 
	end
	ylabel(label)
	title(titlestring)
else 	
	mesh(Tau_1, Tau_2, Data);
	set(gca, 'XScale', xlog, 'YScale', ylog)
	xlabel('\tau_1')
	ylabel('\tau_2')
	zlabel(label)
	axis([min(Tau_1) max(Tau_1) min(Tau_2) max(Tau_2) zmin zmax]);
	title(titlestring)
end
