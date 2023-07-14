function [dist]=binner(ww,wb)
%[dist]=binner(ww,wb) bins 2d data (ww) into bins (wb)


% x=find(ww~=0);
% ntot=length(x);

ntot=size(ww,1)*size(ww,2); %total number of points in population
step=wb(2:end)-wb(1:end-1);

% ww=reshape(ww,[1 ntot]);
% ww=sort(ww);

% wboffset=wb(2:end);
% wb=wb(1:end-1);
% 
% a=find(ww>=wb & ww<wboffset);

for i=2:length(wb)
    %a=find(ww>=wb(i-1) & ww<wb(i));
    %dist(i-1)=length(a)/step(i-1)/ntot; %total number of points in 2d slice in this bin per dw and normalised for all points in slice
    
    a=find(ww>=wb(i-1) & ww<wb(i));
    dist(i-1)=length(a)/step(i-1)/ntot; %total number of points in 2d slice in this bin per dw and normalised for all points in slice
    
end
