function [res,err]=shiftmat(z,n)

err=0;

ndim=ndims(z);

dims=[1:ndim];
dims(find(dims==n))=[];
per=permute(z,[n dims]); %reorder so that the requested dimension is first

siz=size(per);

if siz(1)<3; 
    fprintf(1,'\n***For central differencing need at least 3 elements in differencing dimension***\n');
    err=1;
    break
end

len=siz(1); %length in requested dimension
siz2=siz(2:end);

shift=per(3:end,:);
shift=reshape(shift,[len-2 siz2]);

static=per(1:len-2,:);
static=reshape(static,[len-2 siz2]);

res=shift-static;

% dims2=dims;
% dims2(2:n-1)=dims2(2:n-1)+1;
% 
% if n==1; 
%     dims2=dims;
% else
%     dims2=
% res=permute(z,[n

%need to rearrange so that is in original order