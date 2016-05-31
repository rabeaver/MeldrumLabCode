function mat = orthomult(mat,u,b,Xdim,tailledata)

% global u b Xdim tailledata
% function used to apply the householder transformation to each column in some "mat"
for i = 1 : Xdim
%     c = 0;
    ut=u';
    c = sum((ut(1:tailledata).*mat(1:tailledata,i))./b);
    mat(1:tailledata,i) = mat(1:tailledata,i) + c.*ut(1:tailledata); % to make the algorithm faster
end