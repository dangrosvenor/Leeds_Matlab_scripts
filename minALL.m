function [ma,ii]=minALL(arr)

%finds overall min and indices

si=size(arr);
nodim=length(si);
len=prod(si);

if len==0
    disp('*** No data in the array for minALL !! ***');
end

a=reshape(arr,[1 len]);

[ma i]=min(a);

rem=i;
for i=nodim:-1:2
    block=prod(si(i-1:-1:1)); %size of dimension i blocks
    ii(i)=floor(rem/block); %number of whole blocks contained in len (min=1)
    rem=mod(rem,block); %remainder after whole blocks removed
    if rem~=0
        ii(i)=ii(i)+1; %if is not the end of a row add one
    else
        ii(1:i-1)=si(1:i-1); %if is a multiple of whole block left fill in dimension sizes and finish
        break
    end
end
if rem~=0
    ii(1)=rem; %final amount left is fist index
end
