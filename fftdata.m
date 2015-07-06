function fdata = fftdata(data,willfdata,centers,lengths)
fdata = circshift(data,-(centers-1));
nd = ndims(fdata);
for j=1:nd
    if(willfdata(j))
        fdata = fftshift(fft(fdata,lengths(j),j),j);
    end
end
