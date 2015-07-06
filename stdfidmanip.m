function addata = stdfidmanip(data,dpts,ph0,ap_type,ap_p1,ap_p2,ap_p3)
%Memory hog.  Generates a filter matrix the size of the data set.
%Performs a series of standard manipulations to the data. (no transforms)
% data -> the raw data
%Simple manipulations (Performed first)
% dpts -> drops this number of points from the begining of the FID
% ph0 -> zero order phase adjustment in degrees
%Apodizations.  Each is an array of length = # of dimensions
% ap_type -> 0 = nothing
%               ap_p1,ap_p2,ap_p3 have no effect for these dimensions
%            1 = gaussian
%               apodizes along the given dimension as
%                  exp(-(ap_p1(j-ap_p2)/length)^2)
%            2 = eponential
%               apodizes along the given dimension as
%                   exp(-ap_p1*j/length)

%Get some basic parameters about data.
nd = ndims(data);
shd = size(data);

%Remove the first dpts of each fid.
addata = data(dpts+1:end,:)
shd(1)=shd(1)-dpts
addata = reshape(addata,shd);

%Phase the data (???Am I turning it in the conventional direction???)
addata = addata .* exp(pi*i*ph0/180);

for j=1:nd
    %Make the filter array
    switch ap_type(j)
        case 0
        case 1
            filterarray = exp(-ap_p1(j)*(((1:shd(j))-ap_p2(j))/shd(j)).^2);
        case 2
            filterarray = exp(-ap_p1(j)*(1:shd(j))/shd(j));
    end
    if(ap_type(j) ~= 0)
        %construct the matrix
        if(j > 1)
            filterarray = reshape(filterarray,[ones(1,j-1),shd(j)]);
        else
            filterarray = reshape(filterarray,[shd(j),1]);
        end
        repshd=shd;
        repshd(j)=1;
        bigfilter = repmat(filterarray,repshd);
        
        %apply it
        addata = bigfilter.*addata;
    end
end