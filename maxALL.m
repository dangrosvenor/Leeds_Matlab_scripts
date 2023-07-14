function [ma,ii]=maxALL(arr)
% finds overall max and indices
% [ma,ii]=maxALL(arr)

si=size(arr);
nodim=length(si);
len=prod(si);

a=reshape(arr,[1 len]);

[ma i]=max(a);

rem=i;
for i=nodim:-1:2
    block=prod(si(i-1:-1:1)); %size of dimension i blocks
    ii(i)=floor(rem/block); %number of whole blocks contained in len (min=1)
    rem=mod(rem,block); %remainder after whole blocks removed
    if rem~=0
        ii(i)=ii(i)+1; %if is not the end of a row add one
    else
        ii(1:i-1)=si(1:i-1); %if is a multiple of whole block left fill in dimension sizes
        break
    end
end
if rem~=0
    ii(1)=rem;
end
