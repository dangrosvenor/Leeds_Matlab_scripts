function out=bin2dec_array(bins)
%converts an array of any shape to decimal from binary
%the bits must be in the first dimension of the array



siz=size(bins);
Lsiz=length(siz);


%Expand all the dimensions into a vector form except the first one (the one
%containing the bits)
 bins=bins(:,:);
%this makes the array into an [nbits M] sized array
      
out = squeeze(bins); 

if size(out,2)~=1  %squeezing a [1 1 M] array makes it become of size [M 1] whereas squeezing a [2 1 M]
                   %array gives an array of size [2 M] - so need to
                   %transpose in the latter case
    out=out';
end

out = bin2dec(out);

%reshape to original size (except for the bit dimension) and squeeze
out = squeeze( reshape(out,[siz(2:end)]) );
