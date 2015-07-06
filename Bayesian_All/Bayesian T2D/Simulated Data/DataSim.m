clear
close all
clc

%% Data Simulation Scripts

nrEchoes = 128;
tEcho = 120e-6;   % echo time (s)
T2 = 2e-3; %s?


% T2D Data Sim

gamma = 267.513e6; % Gyromagnetic Ratio (rad T-1 S-1)
G = 6.56; % Magnet Gradient (T/m)
Delta = 2e-3; %Big Delta (s)
D = 1e-08; % Diffusion Coefficient

% Little delta loop
deltaMin = 10e-6; %s
deltaMax = 310e-6;
numSteps = 20; % number of Steps
deltaStep = (deltaMax - deltaMin)/numSteps; % Step Size

data = zeros(nrEchoes, numSteps); % Size of data matrix
count = 0; % step number
%signalNorm = gamma.^2 .* G.^2 .* deltaMin.^2 .* D .* (Delta + (deltaMin./3));
% y=(A).*(exp(-x./t2)); T2 decay
for i = deltaMin:deltaStep:deltaMax
    count = count + 1;
    signal = exp(-gamma.^2 .* G.^2 .* i.^2 .* D .* (Delta + (i./3)));
    for j = 1:nrEchoes
        x = 10;
        while x > 1 || x < -1
            x = randn;
        end
        data(j, count) = (signal .* (exp(-(j.*tEcho)./T2))) + randn/20;
    end
end
dataim = randn(size(data))/20;
data = complex(data, dataim);
data2 = data';

figure(1);
hold on
surf(real(data2)); shading flat;
xlim([1 nrEchoes]); ylim([1 numSteps+1]); zlim([-0.2 1.25]);
xlabel('EchoNum')
ylabel('Delta Number')
zlabel('Signal')
hold off

save('dataSimT2D.dat','data2','-ascii','-tabs')

%% Separate Data by little delta

singleVec = reshape(data, 1, (numSteps+1).*nrEchoes);
firstPntperDelta = zeros(1, numSteps+1);
deltaVector = [deltaMin:deltaStep:deltaMax];

figure(2)
hold on
plot(singleVec);
%xlim([1 5.5e04]);
ylim([-0.20 1.2]);
xlabel('Data Point')
ylabel('Signal')
hold off

figure(3)
hold on
plot(singleVec(1,1:nrEchoes))
xlim([1 nrEchoes]);
ylim([-0.2 1.2]);
xlabel('Echo')
ylabel('Signal')
hold off

for i = 1:(numSteps+1)
    firstPntperDelta(1,i) = data2(i,1);
end

figure(4)
hold on
plot(deltaVector,firstPntperDelta./firstPntperDelta(1,1))
xlim([deltaMin deltaMax])
xlabel('Delta (s)')
ylabel('ln(Signal/S_0)')

fitx = -gamma.^2 .* G.^2 .* deltaVector.^2 .* (Delta + (deltaVector./3));
fity = log(firstPntperDelta./firstPntperDelta(1,1));
figure(5)
hold on
plot(fitx, fity)
ylabel('ln(Signal/S_0)')



