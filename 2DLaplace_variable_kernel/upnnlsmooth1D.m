function [spectrum,tau,chisq,compte]=upnnlsmooth1D(data,time,Tamin,Tamax,alpha,beta,steps,kernel1)

% Non Negative Least Square function with smoothing on the magnitude, slope and curvature of the distribution
% Based on NNLS algorithm of Lawson&Hanson: "solving least squares problems", p 161
% minimization of chisq^2=||Y-BX||^2+alpha||X||^2 + beta||dX/dtau||^2 +beta||d^2X/dtau^2||^2
% where Y is the vector of the noisy data, 
% B the known matrix of exp(-time(i)/Tau(j)) and X the distribution vector 
% statement of the vector tau, depending if the vector data is coming from a 1 or 2 D data set
% If 1D data set, Tbmin = Tbmax = -1

% constrain:
if steps < 4
    warndlg('chose a bigger number of steps');
    return
end

% Definition of the vector tau for the distribution
pas = (log10(Tamax)-log10(Tamin))/(steps-1);
tau = 10.^(log10(Tamin):pas:log10(Tamax));

% constraint
[tailledata,~] = size(data);
[tailletime,~] = size(time);

if tailletime ~= tailledata
    warndlg('error: the time and data vectors must have the same length');
    return
end
pasdata = (log10(time(tailletime,1))-log10(time(1,1)))/(tailletime-1);
% constraint on alpha
constraint = 1/(sqrt(pasdata^2*pas^3*alpha));
if constraint > 300
    warndlg('make alpha bigger');
    return
end

% initialisation for Uniform Penalty loop
time = time';
data = data';
[~,Xdim] = size(tau);

% softcurvature = zeros(1,Xdim);
% softslope = zeros(1,Xdim);
spectrum = zeros(1,Xdim);
beta = beta / (data(1))^2;

if beta < 0 
    Max = 1;
else
    Max = 2;
end

% start of uniform penalty loop
for uploop = 1 : Max
    % definition of the slope and the curvature for the UP loop of the data
    % very first step (or beta = -1) slope = 0 and curvature = 0
    % So only smoothing by alpha (classical one)
    
    slope = diff(spectrum)/pas;
    slope(Xdim) = slope(Xdim-1);
    curvature = diff(diff(spectrum))/pas^2;
    curvature = [0 curvature 0];
    
    softslope = soften(slope,Xdim);
    softcurvature = soften(curvature,Xdim);

    % determination of the weighting which will be used to determine the data
    % first step (or beta = -1): only alpha
    weighting = 1./sqrt(pasdata^2*pas^3*(alpha + beta*(softcurvature.^2 + 5*softslope.^2)));
    
    % number of elements in weighting less than 1e-6
    compte = size(find(weighting < 1e-6));
    
    % constraint
    if Xdim - compte(2) + tailledata - 2 < Xdim
        warndlg('need more data, more curvature constraint or fewer spectrum points');
        return
    end
    
    [spectrum, chisq,compte]=nnlsmooth1Dsvd(spectrum,weighting, tau, time, data, kernel1);
    
    if chisq == -1 
        break
    end
end




