function [s,g,ffit] = reginvlaplace(t,f,m,order,lambda)
% [s,g,ffit] = reginvlaplace(t,f,m,order,lambda)
%
% Salman Rogers (c) 2005 - release version 13 March 2005
%
% Calculates a non-negative, linearly regularised 
% inverse laplace transform of f(t), i.e.:
% f(t)=integral of g(s) with respect to s,
% and plots the results.
%
% Input variables:
%    t and f vector arrays
% Output variables:
%    s and g vector arrays
%    ffit is fitted value of f(t)
% Input parameters:
%    m = number of points s, arranged in a logarithmic scale
%    order = order of linear regularisation: 0, 1, 2 or 3
%    lambda = weight of linear regularisation >0, a lagrange multiplier
%               if lamba=0; no regularisation, 
%               if lambda=infinity; solution is polynomial of order input.
% Recommended: set lambda=-1, then lambda is automatically chosen to 
% make weight of data and regularisation comparable.

n=length(t);
% set output grid s
s_low=1e-2*t(2);
s_high=1e2*t(end);
s=exp(log(s_low) + (log(s_high)-log(s_low))*[0:(1/(m-1)):1]); %log scale
%s=t(2) + (t(end)-t(2))*[0:(1/(m-1)):1]; %linear scale

% calculate K(t,s)=exp(-t/s): matrix a
a=exp(-t*(1./s));

% calculate reg matrix H
ds=diff(s);
switch order
case 0
    B=eye(m); %zeroth order regularisation
case 1
    B=zeros(m-1,m); %first difference matix
    for i=[1:m-1] 
        B(i,i)=-1/ds(i)^0.5;
        B(i,i+1)=1/ds(i)^0.5; %1st order reg
    end
case 2
    B=zeros(m-2,m);
    for i=[1:m-2]
        B(i,i)=-1/ds(i)^1.5;
        B(i,i+1)=1/ds(i)^1.5 + 1/ds(i+1)^1.5;
        B(i,i+2)=-1/ds(i+1)^1.5;
    end
case 3
    B=zeros(m-3,m);
    for i=[1:m-3]
        B(i,i)=-1/ds(i)^2.5;
        B(i,i+1)=1/ds(i)^2.5 + 2/ds(i+1)^2.5;
        B(i,i+2)=-2/ds(i+1)^2.5 - 1/ds(i+2)^2.5;
        B(i,i+3)=1/ds(i+2)^2.5;
    end
end

H=B'*B; 

if lambda=='auto'
    lambda=trace(a'*a)/trace(H)
end

% calculate A = a'.a + lamdba.H
A=a'*a+lambda*H;
% calculate y = a'.b
y=a'*f;
% calculate g by LU decomposition
[L,U]=lu(A);
g=U\(L\y);

% now iterate, mapping negative values to zero.
e=2/max(eig(A));
A=a'*a;  %no regularisation now
for i=[1:1000]
    g=(g>0).*g;    %map neg to zero
    g=(eye(m)-e*A)*g+e*y;
end
g=(g>0).*g;    %map neg to zero again

ffit=a*g;

%convert to differential form
s=(s(1:end-1)+s(2:end))/2;
g=(g(1:end-1)+g(2:end))./ds'/2;

subplot(2,1,1);hold on;plot(t,f);plot(t,ffit,'+r');
xlabel('t'); ylabel('f');
legend('input','fit')
subplot(2,1,2);hold on;plot(s,g,'ob');set(gca,'XScale','log');
xlabel('s'); ylabel('g');
