clear
clc
close all



T2 = 100e-3; %T2 in s
T1 = 0; % T1 in s
nRates = 1; %number of T2 rates
tEcho = 50e-6; %echo time s
tau=tEcho/2; %half the echo spacing
n = 32; %number of echoes
alpha = -180*(2*pi/360);

H1 = zeros(n,nRates);
H2 = zeros(n,nRates);
%
% RF mixing matrix
%
T0 = [cos(alpha/2)^2, sin(alpha/2)^2, sin(alpha); ...
    sin(alpha/2)^2, cos(alpha/2)^2, -sin(alpha); ...
    -0.5*sin(alpha), 0.5*sin(alpha), cos(alpha)];

TArray = cell(1,n);
TArray(:) = {sparse(T0)};
T = blkdiag(1,TArray{:});

%
% Selection matrix to move all traverse states up one coherence level
%
S = sparse(3*n+1,3*n+1);
S(1,3)=1;
S(2,1)=1;
S(3,6)=1;
S(4,4)=1;
for o=2:n
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
  
%

    % Relaxation matrix
    R2 = 1/T2;
    R1 = 1/T1;
    R0 = diag([exp(-tau*R2),exp(-tau*R2),exp(-tau*R1)]);
    
    RArray = cell(1,n);
    RArray(:) = {sparse(R0)};
    R = blkdiag(exp(-tau*R2),RArray{:});

    %Diffusion matrix
    D = 2.2e-9; %m^2 s-1
    G = 23.8626; %Gradient T/m (6.5998 PM25; 23.8626 PM5)
    gamma = 2.675222e8;	% Gamma, Rad/s/T.
    gamma2 = 42.577;    % MHz/T
    k = gamma*G*tau; %k-space traversal in time period tau (tE/2)
    
    DArray = cell(1,n);
    for i = 1:n
        bT = ((i-1)*k+k/2)^2*tau+k^2/12*tau;
        bL = ((i-1)*k)^2*tau;
        Dhold = diag([exp(-D*bT),exp(-D*bT),exp(-D*bL)]);
        DArray{i} = sparse(diag([exp(-D*bT),exp(-D*bT),exp(-D*bL)]));
    end
    
    D2 = blkdiag(1,DArray{:});
    

%RUN EPG INCLUDING DIFFUSION
% Precession and relaxation and diffusion matrix
P = (R*S*D2);

% Matrix representing the inter-echo duration
E = (P*T*P);

% Recursively apply matrix to get echo amplitudes
%
x = zeros(size(R,1),1);
x(1)=1;
    for iEcho=1:n
        x=E*x;
        H2(iEcho,1) = x(1);
    end
    
    
%RUN EPG WITHOUT DIFFUSION
% Precession and relaxation matrix
Pp = (R*S);

% Matrix representing the inter-echo duration
Ee = (Pp*T*Pp);

% Recursively apply matrix to get echo amplitudes
%
x = zeros(size(R,1),1);
x(1)=1;
    for iEcho=1:n
        x=Ee*x;
        H1(iEcho,1) = x(1);
    end


%PLOT COMPARISON OF WITH AND WITHOUT DIFFUSION
figure(1)
hold on
plot(tEcho*(1:n),H1)
plot(tEcho*(1:n),H2)
legend('relaxation only','relaxation diffusion')