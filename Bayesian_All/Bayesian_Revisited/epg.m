% Function: epg
% Author: Jared N. King
% Advisor: Tyler Meldrum
% Source Codes used: Bayesian Estimate Software by Kelvin Layton
% Date Modified: 06/04/2015
% Desc: Program implements an Extended Phase Graph function. First, only
%  the mixing matrix, T, and the selection matrix, S, will be implemented.
%  Second, relaxation matrices, R0 and R1, will be implemented. Third, the
%  feature of Diffusion, noted D, will be implemented.
% Parameters:
%       In: N       - Number of Echoes
%           tau     - echotime
%           R1      - 1/T1
%           R2vec   - 1/T2
%           alpha   - flip angle
%           D       - Diffusion Coefficient
%           G       - Field Gradient (T/m)
%       Out: H      - Representation of signal
%           S       - Selection Matrix (sparse); "moves all transverse 
%               states up one coherence level"
%           T       - RF mixing matrix (block diagonal)
%           R       - Relaxation Matrix (sparse)
%           D       - Diffusion Matrix (sparse)
%           P       - Precession Matrix which will eventually include 
%               relaxation and diffusion(sparse)
%           E       - Matrix representing inter-echo duration (sparse)
%       Returns:    H
% Usage: The function will eventually be used in parallel with a similar
%  derivatives function to implement a Bayesian approach data processer,
%  solving for T2 and Diffusion using data from a T2-D pulse sequence
%  experiement. 
%
%
%   !!! Eventually a function !!!
%   function [H] = epg( N, tau, R1, R2vec, alpha, D, G)
%   !!!                       !!!

clear
close all
clc

% Temp Variables, to be replaced by function call %
n = 4096; % Number of points/echoes
alpha = pi; % Flip angle


    % for use in implementation of T2 in EPG
tau = 100e-6; % echotime
R2vec = [1/(500e-3)]; %#ok<*NBRAK> % 1/T2 (1/s)
R1 = 1/(1500e-3); % 1/T1 (1/s)

% Diffusion Parameters
D = 2.2e-9; % Diffusion m^2 / s
G = 6.56; % Field Gradient (T/m)
gamma = 2.675222e8; % Gamma, Rad/s/T
k = gamma * G * tau; % k-space traversal during time tau

% Diffusion Experiment Extra Params
delta = 10e-06;
DELTA = .05e-03;

%                                                 %

% Function initial variables %
nRates = length(R2vec); % Number of T2 components being solved
tau = tau/2; 

H = zeros(n, nRates); % Signal Matrix Setup; number of echoes by number 
                %of components
    

% RF Mixing Matrix Setup %
%
% Main mixing matrices

% pi pulse
T180 = [cos(pi/2)^2, sin(pi/2)^2, sin(pi); ...
    sin(pi/2)^2, cos(pi/2)^2, -sin(pi); ...
    -0.5*sin(pi), 0.5*sin(pi), cos(pi)];

% pi/2 pulse
T90 = [cos(pi/4)^2, sin(pi/4)^2, sin(pi/2); ...
    sin(pi/4)^2, cos(pi/4)^2, -sin(pi/2); ...
    -0.5*sin(pi/2), 0.5*sin(pi/2), cos(pi/2)];

% Make sparse array the size of n with T0 as diagonal entries
TArray = cell(1,n);
TArray(:) = {sparse(T90)}; % Change main pulse here
TDiff = cell(1, 2);
TDiff(:) = {sparse(T180)}; % For second and third pulse of T2D Experiment

T = blkdiag(1, TDiff{:}, TArray{:}); 


% Selection Matrix Setup %

S = spalloc(3*n + 7, 3*n + 7, 3*n + 6);
% Base
S(1,3) = 1;
S(2,1) = 1;
for i = 3:3*n + 7
    if mod(i,3) == 0
        S(i,i+3) = 1;
    elseif mod(i,3) == 1
        S(i,i) = 1;
    else
        S(i,i-3) = 1;
    end
end

S = S(1:3*n+7, 1:3*n+7);

% Build Signal %

for iRate = 1:nRates % In the case of multiple components
    
    % Relaxation Matrix %
    Rdiff = cell(1,2); % Relaxation during little delta and big DELTA
    RArray = cell(1, n); % pi/2 pulses
    
    R2 = R2vec(iRate); % R2 is 1/T2
    
    % RDELTA is for relaxation during second pulse in sequence
    RDELTA = diag([exp(-DELTA * R2), exp(-DELTA * R2), exp(-DELTA * R1)]);
    
    % Rdelta represents relaxation before 90 pulses begin
    Rdelta = diag([exp(-delta * R2), exp(-delta * R2), exp(-delta * R1)]);
    
    % During 90 pulses
    R0 = diag([exp(-tau * R2), exp(-tau * R2), exp(-tau * R1)]); % Matrix 
                                                    % for block diagonal
    % Create Relaxation matrix
    RArray(:) = {sparse(R0)}; 
    R = blkdiag(exp(-tau * R2), RDELTA, Rdelta, RArray{:}); % This is R
                                                    
    % Diffusion Matrix %
    
    % first pulse = 1
    % Second pulse, big delta
    % Third pulse, little delta
    Ddiff = cell(1, 2); % Matrix for second and third pulses
    DArray = cell(1,n);
    
    for i = 1:n + 2 % n for nechoes and 2 for extra 90 pulses
        if i == 1
            %Transverse Relaxation during big DELTA
            bT = ((i-1) * k + k/2)^2 * DELTA + k^2/12 * DELTA;
            
            % Longitudinal Relaxation during big DELTA
            bL = ((i-1) * k)^2 * DELTA; 
            
            % Add to matrix
            Ddiff{1} = sparse(diag([exp(-D * bT), exp(-D*bT), exp(-D*bL)]));
        
        
        elseif i == 2
            %Transverse Relaxation during little delta
            bT = ((i-1) * k + k/2)^2 * delta + k^2/12 * delta;
            
            % Longitudinal Relaxation during little delta
            bL = ((i-1) * k)^2 * delta;
            
            % Add to matrix
            Ddiff{2} = sparse(diag([exp(-D * bT), exp(-D*bT), exp(-D*bL)]));
            
        else
        
            bT = ((i-1) * k + k/2)^2 * tau + k^2/12 * tau; %Transverse Relaxation
                                                        % at each point
            bL = ((i-1) * k)^2 * tau; % Longitudinal Relaxation at each point
            
            % Sparse Diagonal of effective diffusion at each point
            DArray{i-2} = sparse(diag([exp(-D*bT), exp(-D*bT), exp(-D*bL)]));
        end
    end
    
    
        
    DMatrix = blkdiag(1, Ddiff{:}, DArray{:}); % Block Diagonal Diffusion Matrix
                                % [Default, 90 pulse, 90 pulse, 180...n]
                
    
    % Precession %
    P = (R * DMatrix * S);
           
    % Inter-Echo Duration Matrix %
    E = (P*T*P); % Pulse, echotime, pulse
    
    % Recursively apply E to get echo amplitudes %
    x = zeros(size(S, 1),1);
    x(1) = 1;
    for iEcho = 1:n
        x = E*x;
        H(iEcho, iRate) = x(1);
    end
end

% end %

