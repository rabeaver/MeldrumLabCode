function [setp,answ1]=lineareqsolve(systemmat,setp,answ1,Xdim1)

% systemmat has already been diagonalized, only need to know Xdim

ind = find(setp == 0);
setlength = ind(1)-1;

tempans = zeros(1,setlength);       % initialise tempans
tempans(setlength) = systemmat(setlength,Xdim1)/systemmat(setlength,setp(setlength));
for i = setlength-1 :-1 : 1
%     som = 0;
    som = sum(systemmat(i,setp(setlength:-1:i+1)).*tempans(setlength:-1:i+1));
    tempans(i) = (systemmat(i,Xdim1)-som)/systemmat(i,setp(i));    
end

answ1=zeros(1,Xdim1-1);  
answ1(setp(1:setlength))=tempans(1:setlength);
