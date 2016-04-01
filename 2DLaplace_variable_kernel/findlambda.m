function [lambda,zeroindice]=findlambda(spectrum,soln,setp,Xdim)

% step #8
% Find an indice q (zeroindice) in setp such that xq/(xq-zq)=min
% x=spectrum and z = soln (solution of step 6)
zeroindice=0;
lambda=0;
for i = 1 : Xdim
    if setp(i) ~= 0
%        sol = soln(setp(i));
        
       if soln(setp(i))<0
            % run thought all indices in setp, 
            %calculate the lambda value for indices that correspond 
            %to negative coefficients in soln, see which lambda is smallest
            if lambda == 0
                lambda = spectrum(setp(i))/(spectrum(setp(i))-soln(setp(i)));
                % for the first time in the loop
                zeroindice = setp(i);
            end
            if spectrum(setp(i))/(spectrum(setp(i))-soln(setp(i))) < lambda
%                 la = 2;
                lambda = spectrum(setp(i))/(spectrum(setp(i))-soln(setp(i)));
                % successive loops compared to current lambda
                zeroindice = setp(i);
            end
        end
    end
end
% step # 9: set lambda = xq/(xq-zq)