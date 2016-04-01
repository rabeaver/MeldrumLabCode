function [systemmat] = qrdecmp(systemmat,setp,oldsetp,tailleX,tailleY)

% series of Householder tranformation designed sepcifically for the LS problem

setlength = 1;
added = 1;
removed = 1;
u=zeros(1,tailleY);

while setp(setlength)~=0 && oldsetp(setlength)~=0
    setlength = setlength + 1;
end
% find a zero to either setp or oldsetp

if setp(setlength)==0
    added = 0;
end

if added == 1
    % householder transformation to most recently chosen column, setp(setlenght) in the matrix
    [b,u] = house(setp(setlength),setlength,systemmat,tailleY,u);
    systemmat = orthomult(systemmat,u,b,tailleX,tailleY);
else
    % householder transformation to all columns listed in setp
    while setp(removed) == oldsetp(removed)
        removed = removed + 1;
    end
    k = removed;
    while setp(k) ~= 0
        [b,u] = house(setp(k),k,systemmat,tailleY,u);
        systemmat = orthomult(systemmat,u,b,tailleX,tailleY);
        k = k+1;
    end
end
