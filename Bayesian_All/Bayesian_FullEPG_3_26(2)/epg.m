% EPG Matrix implementation of the Extended Phase Graph (EPG) algortihm.
%   H = EPG(N,TAU,R1,R2,ALPHA) calculate the echo amplitudes for a CPMG
%   sequence with N echoes and echo spacing of (2*TAU). The parameters are
%   R1=1/T1 (scalar) and R2=1./T2 (scalar or vector) and ALPHA the flip 
%   angle. The matrix H contains one column for each element in R2
%
% Author: Kelvin Layton (klayton@unimelb.edu.au)
% Date: June 2012
%
function [ H ] = epg( N, tau, R1, R2vec, alpha , D, G)
nRates = length(R2vec);
tau=tau/2;

H = zeros(N,nRates);

% Bvals and Diffusion thing
% bvalm = 1.0428e6;
% bvalp = 1.0428e6;
% bvalZ = 0;
% D = 1.0e-11;

% RF mixing matrix
%
T0 = [cos(alpha/2)^2, sin(alpha/2)^2, sin(alpha); ...
    sin(alpha/2)^2, cos(alpha/2)^2, -sin(alpha); ...
    -0.5*sin(alpha), 0.5*sin(alpha), cos(alpha)];

TArray = cell(1,N);
TArray(:) = {sparse(T0)};
T = blkdiag(1,TArray{:});

% Selection matrix to move all traverse states up one coherence level
%
S = sparse(3*N+1,3*N+1);
S(1,3)=1;
S(2,1)=1;
S(3,6)=1;
S(4,4)=1;
for o=2:N
    offset1=( (o-1) - 1)*3 + 2;
    offset2=( (o+1) - 1)*3 + 3;
    if offset1<=(3*N+1)
    S(3*o-1,offset1)=1;  % F_k <- F_{k-1}
    end
    if offset2<=(3*N+1)
    S(3*o,offset2)=1;  % F_-k <- F_{-k-1}
    end
    S(3*o+1,3*o+1)=1;  % Z_order
end


% gamma = 2.675222e8;	% Gamma, Rad/s/T.
% % gamma2 = 42.577;    % MHz/T
% 
% 
% k  = gamma*G*tau;
% D0 = diag([exp(-7*k^2*tau/3*D), exp(-7*k^2*tau/3*D), exp(-k^2*tau*D)]);
% % D0d = diag([-7*k^2*tau/3*exp(-7*k^2*tau/3*D), -7*k^2*tau/3*exp(-7*k^2*tau/3*D), -k^2*tau*exp(-k^2*tau*D)]);
% DArray = cell(1,N);
% DArray(:) = {sparse(D0)};
% D_mat = blkdiag(exp(-7*k^2*tau/3*D),DArray{:});


for iRate=1:nRates

    % Relaxation matrix
    R2=R2vec(iRate);
    R0 = diag([exp(-tau*R2),exp(-tau*R2),exp(-tau*R1)]);

    RArray = cell(1,N);
    RArray(:) = {sparse(R0)};
    R = blkdiag(exp(-tau*R2),RArray{:});
    
     % Diffusion Matrix
%     D0 = diag([exp(-bvalp*D),exp(-bvalm*D),exp(-bvalZ*D)]);
%     
%     DArray = cell(1,N);
%     DArray(:) = {sparse(D0)};
%     D1 = blkdiag((bvalp*D),DArray{:});

    % Precession and relaxation matrix
%     P = (D_mat*R*S);
    P = (R*S);

    % Matrix representing the inter-echo duration
    E = (P*T*P);
    
    % Recursively apply matrix to get echo amplitudes
    %
    x = zeros(size(R,1),1);
    x(1)=1;
    for iEcho=1:N
        x=E*x;
        H(iEcho,iRate) = x(1);
    end

end
end

