%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                            %
%    Minimize function f Called by Estimate.                                 %
%    The function outputs the function value, the gradient and the Hessian   %
%                                                                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [f,g, H] = minfun(c, data,  K, alpha,Identity)

   % Compute matrix M(c)
   M = ComputeM(K,c);
   % Compute M2 = M + alpha*Identity 
   M2 = alpha*Identity + M;
   % Temporary matrix storage 
   M3 = M2*c;
   % Function to be minimzed
   f = 0.5 * c' * M3 - c'*data;
   % Gradient of the function
   if (nargout > 1) g = M3-data; end
   % Hessian of the function
   if (nargout > 2) H = M2; end



return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                  %
%   Compute matrix G(c), called during calculation of function     %
%   to be minimized. This is a crucial time-consuming step.        %
%                                                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [M] = ComputeM(K, C)

  V = (K'*C > 0); G = K;
  % Construct vector i whose elements are indices of V whose elements are zero
  i = find(V == 0); G(:,i) = 0;
  M = G*K';

return
