function [setp,setz]=deletefromset(setp,setz,index,Xdim)

% removes element in setp that has spectrum(setp) = 0
setz(index) = 1;
i = 1;
while setp(i) ~= index
    i=i+1;
end
for j = i : Xdim-1
    setp(j) = setp(j+1);
end
setp(Xdim)=0;
