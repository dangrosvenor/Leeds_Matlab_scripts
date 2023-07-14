function [x,y] = window_average(xdat,ydat,Nwindow,operation_flag)
%function [x,y] = window_average(xdat,ydat,Nwindow,operation_flag)
%For an Nd array it assumes that the FIRST dimension is the time dimension
%- others are collapsed and restored after.
siz=size(ydat);
nd=length(find(siz>1)); %Number of non-singleton dimensions
if nd>1
    ydat=ydat(:,:); %collapse into 2 dim array
end


switch operation_flag
    case 'mean'
        bfilter=ones([1 Nwindow])*1/Nwindow;
    case 'sum'
        bfilter=ones([1 Nwindow])*1;
end

%do the window operation - this will give an array the same size as y -
%however, the first 1:Nwindow-1 entries will not use the full Nwindow
%datapoints and aren't required for our purposes.
y=filter(bfilter,1,ydat);

%remove the entries not requried

if nd==1
    y(1:Nwindow-1)=[];
else
    y(1:Nwindow-1,:)=[];
    y = reshape(y,[size(y,1) siz(2:end)]);
end


%now fix x so that it corresponds to the mid-points of the data over which
%the window was applied
if mod(Nwindow,2)==0  %even number - in this case the window mid-points falls between datapoints in x
    %and so we will find all of the points halfway between 
    x =0.5*(xdat(Nwindow/2:end-Nwindow/2) + xdat(Nwindow/2+1:end-Nwindow/2+1));
else %for odd Nwindow values we can just use the exact x datapoints starting at the mid-point of the first window
    x = xdat(Nwindow/2+0.5:end-Nwindow/2+0.5);
end

