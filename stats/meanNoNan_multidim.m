function [me,nnums]=meanNoNan_multidim(dat,ndim)
%function [me,nnums]=meanNoNan_multidim(dat,ndim)
%finds the mean of an N dimensional array without the NaNs over dimension ndim
%nnums is the number of values that went into each mean value (no. non-NaN
%values)


%dat=squeeze(dat);
inotnan=~isnan(dat); %so returns a matrix of ones for points that aren't NaN

size_dat=size(dat);

% sum_dims=1:length(size_dat); %list of all the dimensions of dat
% sum_dims(sum_dims==ndim)=[]; %remove the dimension over which we require the mean
% 
% sum_dat=dat;
% for i=1:length(sum_dims)
%     sum_dat=sum(sum_dat,sum_dims(i));  %sum over the other dimensions
%     
    
nnums=sum(inotnan,ndim); %number of non-NaN numbers in the ndim dimension
ib=find(inotnan==0); %all the points that are NaN
dat(ib)=0; %set them to zero


me=sum(dat,ndim)./nnums; %multiplication by inan removes all NaN numbers from sum

me(nnums==0)=NaN; %set the answer to NaN when have no non-NaN data for an entry

me=squeeze(me);