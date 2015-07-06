%This routine implements Hankel transforms of integer order based on a
%Fourier-Bessel series expansion. The algorithm is based on a recently published 
%research work (please cite if used):

%M. Guizar-Sicairos and J. C. Gutierrez-Vega, Computation of quasi-discrete
%Hankel transforms of integer order for propagating optical wave fields,
%J. Opt. Soc. Am. A 21, 53-58 (2004).

%The numerical method features great accuracy and is energy preserving by
%construction, it is especially suitable for iterative transformation
%processes. Its implementation, requires the computation of zeros of 
%the m-th order Bessel function of the first kind where m is the
%transformation order.

%An array of the first 3001 Bessel functions of order from cero to four
%can be found in the "c.mat" array. If a greater transformation order is
%required the zeros may be found numerically. 
%With the c.mat array, as included, Hankel transforms of order 0-4 may be
%computed, with up to 3000 sampling points. A trasformation, and inverse 
%transformation example is given below.
% This routine was tested under Matlab 6.5 R13

% clear
close all
clc

%
%%%%%%%%%%%%%%%%%%%%%%%
%%Input parameters  %%
%%%%%%%%%%%%%%%%%%%%%%%
R = max(profile1(:,1))
N = size(profile1,1)
% R = 150          %% Maximum sampled radius (time)
% N = 1000        %% Number of sampling points
ord = 1        %% Transformation order

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Matrix and vectors computing  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This operations may only be performed once for iterative algorithms

load C:\Users\TMeldrum\Documents\MATLAB\c.mat; %load in a matrix with the zeroes fo the Bessel functions
c = c(ord+1,1:N+1); %set a vector 'c' of the zeroes of the next (n+1) Bessel for each sampled point

V = c(N+1)/(2*pi*R);    % Maximum frequency
r = profile1(:,1);
% r = c(1:N)'*R/c(N+1);   % Radius vector
v = c(1:N)'/(2*pi*R);   % Frequency vector

[Jn,Jm] = meshgrid(c(1:N),c(1:N)); %make two arrays (N x N), one with the zeroes in every row and one with them in every column
C = (2/c(N+1))*besselj(ord,Jn.*Jm/c(N+1))./(abs(besselj(ord+1,Jn)).*abs(besselj(ord+1,Jm)));
% C is the transformation matrix

m1 = (abs(besselj(ord+1,c(1:N)))/R)';   %% m1 prepares input vector for transformation
m2 = m1*R/V;                            %% m2 prepares output vector for display
clear Jn
clear Jm

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Transformation example  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The n-th order Hankel transform of a generalized top-hat function is obtained and
%compared with the analytic solution, please chage transformation order
%from 0 to 4 to see additional results.

f = profile2(:,1);
% f = (r.^ord).*(r<=1);        %% Generalized top-hat function (input)
F = f./m1;                   %% Prepare vector for transformation
F2 = C*F;                    %% Obtain the Hankel transform
Fretrieved = C*F2;           %% Inverse hankel transform

fretrieved = Fretrieved.*m1; %% Prepare vector for display
f2 = F2.*m2;                 %% Prepare vector for display

%%
%%%%%%%%%%%%%%%%%%%%%%%%
%Display results  %%
%%%%%%%%%%%%%%%%%%%%%%%%

figure(1), 
subplot(3,2,1), plot(r,f), 
xlabel('r'),
title('input function')
subplot(3,2,2), plot(r,fretrieved,'r'),
xlabel('r'),
title('Retrieved function with IHT'),
subplot(3,2,3:4), plot(r,f - fretrieved)
xlabel('r')
legend('residuals')
v2 = linspace(1e-10,5,300);
fanalytic = besselj(ord+1,2*pi*v2)./v2;
subplot(3,2,5:6), plot(v,f2,'.r'), hold on, plot(v2,fanalytic), hold off, xlim([0 5]),
legend('Transformation results','Analytic Solution'),
xlabel('v')

% All codes were written by Manuel Guizar Sicairos.