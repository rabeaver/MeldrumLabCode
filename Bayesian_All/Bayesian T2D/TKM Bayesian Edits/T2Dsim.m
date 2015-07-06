function [simSpec] = T2Dsim(nrEchoes, tEcho, T2, G, DELTA, D, deltaMin, deltaMax, n2DPts, noiseFactor)

gamma = 267.513e6; % Gyromagnetic Ratio (rad T-1 S-1)
deltaVec = linspace(deltaMin,deltaMax,n2DPts); % vector of delta times (s)

data = zeros(nrEchoes, n2DPts); % Size of data matrix
for i = 1:n2DPts
    signal = exp(-gamma.^2 .* G.^2 .* deltaVec(i).^2 .* D .* (DELTA + (deltaVec(i)./3)));
    for j = 1:nrEchoes
        data(j, i) = (signal .* (exp(-(j.*tEcho)./T2))) + randn/noiseFactor;
    end
end

simSpec = data';
end