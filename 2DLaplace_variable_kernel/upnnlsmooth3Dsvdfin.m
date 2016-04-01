function [spectrum,taua,taub,chisq,compte]=upnnlsmooth3Dsvdfin(data,timea,timeb,Tamm,stepsa,Tbmm,stepsb,alpha,beta,orient,kernel1,kernel2)

% Non Negative Least Square function with smoothing on the magnitude, slope and curvature of the distribution
% Based on NNLS algorithm of Lawson&Hanson: "solving least squares problems", p 161

% minimization of chisq^2=||Y-BX||^2+alpha||d^2X/dtau^2||^2 (+ beta||dX/dtau||^2 +beta||d^2X/dtau^2||^2)
% where Y is the vector of the noisy data, 
% B the known matrix of exp(-time(i)/Tau(j)) and X the distribution vector 
% last modification: adaptation of the algorithm to a 2D data set:
% minimization of ||Y-B1XB2'||^2 + alpha ||X||^2
% First step : determination of B1 and B2, svd applied on each matrice (B1, B2, data)
% and then and only then: reshape of the data matrix to a 1D on
% and determination of the B matrix product of B1 and B2


% statement of the vector tau, depending if the vector data is coming from a 1 or 2 D data set
% If 1D data set, Tbmin = Tbmax = -1

% constrain:
if stepsa < 4 || (stepsb < 4 && stepsb ~= 1)
    warndlg('chose a bigger number of steps');
    return
end

% Definition of the vector tau for the distribution
Tamin = Tamm(1);
Tamax = Tamm(2);
Tbmin = Tbmm(1);
Tbmax = Tbmm(2);
% taua = zeros(1,stepsa);
pasa = (log10(Tamax)-log10(Tamin))/(stepsa-1);
taua = 10.^(log10(Tamin):pasa:log10(Tamax));

if stepsb ~= 1
%     taub = zeros(1,stepsb);
    pasb = (log10(Tbmax)-log10(Tbmin))/(stepsb-1);
    taub = 10.^(log10(Tbmin):pasb:log10(Tbmax));
else
    taub = 1;
end

% constraint to avoid wrong data and time set
[tailledatab,tailledataa] = size(data);
ttimea = length(timea);
ttimeb = length(timeb);
if ttimea ~= tailledataa && ttimeb ~= tailledatab
    warndlg('error: the times and data matrices must have the same dimension');
    return
end



% constraint on alpha
pas = pasa; % Constraint applied on the dimension "a"
pasdata = (log10(timea(ttimea))-log10(timea(1,1)))/(ttimea-1);
constraint = 1/(sqrt(pasdata^2*pas^3*alpha));
if constraint > 300
    warndlg('make alpha bigger');
    return
end

% initialisation for Uniform Penalty loop

[~,Xadim] = size(taua);
[~,Xbdim] = size(taub);

% Xdim = Xadim*Xbdim;
% softcurvature = zeros(1,Xadim*Xbdim);
% softslope = zeros(1,Xadim*Xbdim);
spectrum = zeros(1,Xadim*Xbdim);
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
    % So only smoothing by alpha (classical one), on the curvature
    
    slope = diff(spectrum)/pasa;
    slope(Xadim*Xbdim) = slope(Xadim*Xbdim-1);
    curvature = diff(diff(spectrum))/pasa^2;
    curvature = [0 curvature 0];
    
    softslope = soften(slope,Xadim*Xbdim);
    softcurvature = soften(curvature,Xadim*Xbdim);

    % determination of the weighting which will be used to determine the data
    % first step (or beta = -1): only alpha
    weighting = 1./sqrt(pasdata^2*pasa^3*(alpha + beta*(softcurvature.^2 + 5*softslope.^2)));
    
    % number of elements in weighting less than 1e-6
    compte = size(find(weighting < 1e-6));
    
%     new size of the data after compression of the data but svd
    tailledata = tailledataa*tailledatab;
%     Xdim=Xadim*Xbdim
%     co=compte(2)
    % constraint on the number of weighting 
    if Xadim*Xbdim - compte(2) + tailledata - 2 < Xadim*Xbdim
        warndlg('need more data, more curvature constraint or fewer spectrum points');
        return
    end
    
    if orient == 'h'
        [spectrum, chisq,compte]=nnlsmooth3Dsvdhorfin(data, timea, timeb, taua, taub, spectrum, weighting, kernel1, kernel2);
    elseif orient == 'v'
        [spectrum, chisq,compte]=nnlsmooth3Dsvdverfin(data, timea, timeb, taua, taub, spectrum, weighting, kernel1, kernel2);
    elseif orient == 'b'
        [spectrumh, chisq,compte]=nnlsmooth3Dsvdhorfin(data, timea, timeb, taua, taub, spectrum, weighting, kernel1, kernel2);
		compte_temp=compte;
        [spectrumv, chisq,compte]=nnlsmooth3Dsvdverfin(data, timea, timeb, taua, taub, spectrum, weighting, kernel1, kernel2);
		compte=compte+compte_temp;
    end
    
    if chisq == -1 
        break
    end
end
% resize spectrum from a 1D data set to the 2D data set
if orient == 'h'
    spectrum = reshape(spectrum,Xadim,Xbdim);
    spectrum = spectrum';
elseif orient == 'v'
    spectrum = reshape(spectrum,Xbdim,Xadim);
elseif orient == 'b'
    spectrumh = reshape(spectrumh,Xadim,Xbdim);
    spectrumh = spectrumh';
    spectrumv = reshape(spectrumv,Xbdim,Xadim);
    spectrum = sqrt(spectrumh .* spectrumv);
end
