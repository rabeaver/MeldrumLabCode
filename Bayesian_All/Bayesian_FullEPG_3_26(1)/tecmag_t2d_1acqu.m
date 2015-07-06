

% Tecmag-to-matlab T2D Data processing
% Output give [num2Dpnts x nrEchoes]
function [aq, spec2] = tecmag_t2d_1acqu(datafile, nrEchoes, tEcho, pw, d4)
[aq, ~, spec2] = readTecmag4d(datafile);


catArray = zeros(length(spec2(:,1)), pw + d4); % number of zeroes to add to the end
spec2 = cat(2, spec2, catArray);
spec2 = reshape(spec2, length(spec2(:,1)), tEcho, nrEchoes);
spec2 = sum(spec2,2);
spec2 = reshape(spec2, length(spec2(:,1)), nrEchoes);
surf(real(abs(spec2))); shading flat