
function  [FEst, CompressedData, Chi, alpha] = FLIEstimate(Data, Tau_1, ...
Tau_2,  U1, U2, V1, V2,S1, S2,  T1, T2, alpha, NoiseStd, flag, scale)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                %
%   Function is called by the main function FLI.m                                             %
%   Estimate of F(T1, T2), subject to F(T1, T2) >=0, from the data. 
%                                                                                                %
%   Function inputs are :                                                                        %
%   1. Data, details about the data                                                              %
%   2. SVD of matrices K_2 and K_1                                                               %
%   3. Values of alpha for which we need to solve the problem                                    %
%   4. Flag = 1 : if alpha automaticaly chosen (alphas = alpha_start)                            %
%      Flag = 0 : if alpha fixed at alphas = alpha_fixed                                         %
%   Function outputs are :                                                                       %
%   1. FEst for a given value of $\alpha$.                                                       %
%   2. Compressed and Projected data                                                             %
%   3. Optimum value of \alpha (if flag == 0)                                                    %
%                                                                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
