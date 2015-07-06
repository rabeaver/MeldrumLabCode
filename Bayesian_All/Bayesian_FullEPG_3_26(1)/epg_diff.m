% EPG Matrix implementation of the Extended Phase Graph (EPG) algortihm.
%   H = EPG(N,TAU,R1,R2,ALPHA) calculate the echo amplitudes for a CPMG
%   sequence with N echoes and echo spacing of (2*TAU). The parameters are
%   R1=1/T1 (scalar) and R2=1./T2 (scalar or vector) and ALPHA the flip 
%   angle. The matrix H contains one column for each element in R2
%
% Author: Kelvin Layton (klayton@unimelb.edu.au)
% Date: June 2012
%
% Modified for SSE diffusion, T Meldrum, 2015

function [ H ] = epg_diff( n, tau, R1, R2vec, alpha, D, G, delta, DELTA)

nRates = length(R2vec);
tau=tau/2;

H = zeros(n,nRates); %added n+2 to acommodate SSE pulse seq

% RF mixing matrix
%
T0 = [cos(alpha/2)^2, sin(alpha/2)^2, sin(alpha); ...
    sin(alpha/2)^2, cos(alpha/2)^2, -sin(alpha); ...
    -0.5*sin(alpha), 0.5*sin(alpha), cos(alpha)];
T1 =  [cos(alpha/4)^2, sin(alpha/4)^2, sin(alpha/2); ...
    sin(alpha/4)^2, cos(alpha/4)^2, -sin(alpha/2); ...
    -0.5*sin(alpha/2), 0.5*sin(alpha/2), cos(alpha/2)];

% Build Matrix for All States
TArray = cell(1,n);
TArray(:) = {sparse(T0)};
TArray1 = cell(1,2);
TArray1(:) = {sparse(T1)};
T = blkdiag(1,TArray1{:},TArray{:});

% Selection matrix to move all traverse states up one coherence level
%
S = sparse(3*n+1,3*n+1);
S(1,3)=1;
S(2,1)=1;
S(3,6)=1;
S(4,4)=1;
for o=2:n+2
    offset1=( (o-1) - 1)*3 + 2;
    offset2=( (o+1) - 1)*3 + 3;
    if offset1<=(3*n+1)
    S(3*o-1,offset1)=1;  % F_k <- F_{k-1}
    end
    if offset2<=(3*n+1)
    S(3*o,offset2)=1;  % F_-k <- F_{-k-1}
    end
    S(3*o+1,3*o+1)=1;  % Z_order
end
    
for iRate=1:nRates

    % Relaxation matrix
    R2=R2vec(iRate);
    R0 = diag([exp(-tau*R2),exp(-tau*R2),exp(-tau*R1)]);
    Rdelta = diag([exp(-delta*R2),exp(-delta*R2),exp(-delta*R1)]); %relaxation during delta, DELTA periods
    RDELTA = diag([exp(-DELTA*R2),exp(-DELTA*R2),exp(-DELTA*R1)]);

    RArray = cell(1,n);
    RArray(:) = {sparse(R0)};
    R = blkdiag(exp(-delta*R2),RDELTA,Rdelta,RArray{:}); %added exp(-delta...) and RDELTA, Rdelta to accommodate pulse seq

    % Diffusion matrix
    gamma = 2.675222e8;	% Gamma, Rad/s/T.
    k = gamma*G*tau; %k-space traversal in time period tau (tE/2)
    kdelta = gamma*G*delta; %pulse seq mod
    kDELTA = gamma*G*DELTA; %pulse seq mod
%     
    bTDELTA         = ((1-1)*kDELTA+kDELTA/2)^2*DELTA+kDELTA^2/12*DELTA; 
    bLDELTA         = ((1-1)*kDELTA)^2*DELTA; 
    DArrayDELTA{:}  = sparse(diag([exp(-D*bTDELTA),exp(-D*bTDELTA),exp(-D*bLDELTA)]));

% Need to look up more about how to use bT, bL and diffusion for EPG matrix
% calculations. See Weigel.

%    DArrayDELTA{:}  = sparse(diag([1,1,1]));
    
    bTdelta         = ((2-1)*kdelta+kdelta/2)^2*delta+kdelta^2/12*delta; 
    bLdelta         = ((2-1)*kdelta)^2*delta; 
    DArraydelta{:}  = sparse(diag([exp(-D*bTdelta),exp(-D*bTdelta),exp(-D*bLdelta)]));
    
    DArray = cell(1,n);
    for i = 3:n+2
        bT = ((i-1)*k+k/2)^2*tau+k^2/12*tau; 
        bL = ((i-1)*k)^2*tau; 
        DArray{i-2} = sparse(diag([exp(-D*bT),exp(-D*bT),exp(-D*bL)]));
    end
    
    D2 = blkdiag(1,DArrayDELTA{:},DArraydelta{:},DArray{:});
    
    % Precession and relaxation matrix
    P = (R*D2*S);
%     P2 = (R*S*D2);
%     P3 = (S*R*D2);
%     P4 = (S*D2*R);
%     P5 = (D2*S*R);
%     P6 = (D2*R*S);    

    % Matrix representing the inter-echo duration
    E = (P*T*P);
    
    % Recursively apply matrix to get echo amplitudes
    %
    x = zeros(size(R,1),1);
    x(1)=1;
    for iEcho=1:n+2
        x=E*x;
        if iEcho >= 3
        H(iEcho-2,iRate) = x(1);
        end
%         H(iEcho,iRate,2) = x(2);
%         H(iEcho,iRate,3) = x(3);
%         H(iEcho,iRate,4) = x(4);
%         H(iEcho,iRate,5) = x(5);
%         H(iEcho,iRate,6) = x(6);
%         H(iEcho,iRate,7) = x(7);
%         H(iEcho,iRate,8) = x(8);
%         H(iEcho,iRate,9) = x(9);
%         H(iEcho,iRate,10) = x(10);
    end

end

end

