function [me,nnums,stdev]=meanNoNan(dat,ndim,op,isqueeze,iweight,weights)
%function [me,nnums,stdev]=meanNoNan(dat,ndim,op,isqueeze,iweight,weights)
%finds the mean of an N dimensional array without the NaNs over dimension ndim
%nnums is the number of values that went into each mean value (num of non-NaN
%values). And std dev. For sum set op='sum'
%If iweight is set then it will do a weighted mean using the weights given
% I.e weighted_mean = sum( x_i * weight_i) / sum(weight_i)
% Set op to 'sum' to do the sum instead of the mean

if nargin<3
    op='';
end

if nargin<4
    isqueeze=1;
end

if nargin<5
    iweight=0;
    weights = ones(size(dat));
end

sdat=size(dat);

%Special case where data with a zero dimension has been input
i0=find(sdat==0);
if length(i0)>0
    sdat2=sdat;
    sdat2(i0)=[];
    me = NaN*ones(sdat2);
    nnums = me;
    stdev = me;
    fprintf(1,'*** WARNING - zero dimension in dat in meanNaNan - setting output to NaNs');
    return
end
      
    
dims=1:length(sdat);
if length(dims)<ndim
    me=dat;
    nnums=1;
    stdev=NaN;
    return
else
    dims(ndim)='';
end

%Put the averaging dimension (ndim) first.
dat = permute(dat,[ndim dims]);
weights = permute(weights,[ndim dims]);

%convert to a 2D array of single vectors of the length of the ndim
%dimension
dat = dat(:,:);
weights = weights(:,:);

%Chunk data into sections if very large
ntot = prod(size(dat(1,:)));
max_chunk = 5e8;
[start_inds,end_inds]=chunk_inds(ntot,max_chunk);

for i=1:length(start_inds)
    
    [me_chunk,nnums_chunk,stdev_chunk]=meanNoNan_FUNC(dat(:,start_inds(i):end_inds(i)),1,op,0,iweight,weights(:,start_inds(i):end_inds(i)));
    
    if i==1
        me = me_chunk;
        nnums = nnums_chunk;
        stdev = stdev_chunk;        
    else
        me = cat(1,me,me_chunk);
        nnums = cat(1,nnums,nnums_chunk);
        stdev = cat(1,stdev,stdev_chunk);
    end
    
end

sdat_new = sdat;
sdat_new(ndim)=[];
if length(sdat_new)>1
    me = reshape(me,sdat_new);
    nnums = reshape(nnums,sdat_new);
    stdev = reshape(stdev,sdat_new);
end


 if isqueeze==1
    me=squeeze(me);  
    nnums=squeeze(nnums);
    stdev=squeeze(stdev);
 end
    




function [me,nnums,stdev]=meanNoNan_FUNC(dat,ndim,op,isqueeze,iweight,weights)
%function [me,nnums,stdev]=meanNoNan(dat,ndim,op,isqueeze,iweight,weights)
%finds the mean of an N dimensional array without the NaNs over dimension ndim
%nnums is the number of values that went into each mean value (num of non-NaN
%values). And std dev. For sum set op='sum'
%If iweight is set then it will do a weighted mean using the weights given
% I.e weighted_mean = sum( x_i * weight_i) / sum(weight_i)
% Set op to 'sum' to do the sum instead of the mean

if nargin<4
    isqueeze=1;
end

if nargin<5
    iweight=0;
    weights = ones(size(dat));
end

sdat=size(dat);
dims=1:length(sdat);
if length(dims)<ndim
    me=dat;
    nnums=1;
    stdev=NaN;
    return
else
    dims(ndim)='';
end

%Put the averaging dimension (ndim) last.
dat = permute(dat,[dims ndim]);
weights = permute(weights,[dims ndim]);

ndim=length(sdat);


%dat=squeeze(dat);
inotnan=~isnan(dat); %returns a matrix of ones for points that aren't NaN

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
%clear inotnan
dat(ib)=0; %set them to zero
weights(ib)=0;

if iweight==0
    me=sum(dat,ndim)./nnums; %
else
   nweights = sum(weights,ndim);
   me=sum(dat.*weights,ndim)./nweights;  
end

n_ones = ones([1 length(sdat)-1]);

me2 = repmat(me,[n_ones size(dat,ndim)]);

if iweight==1
    %In this case want sum( (w_i * (x_i - weighted_mean).^2) ) /
    %                           sum( N/N-1 * w_i)
    % where N is the number of non-zero weights (NIST formula)
    
    %weights has been set to zero where there is NaN data in dat
    A =  sum( weights .* (dat - me2).^2 , ndim);
    B = (nweights-1)./(nweights) .* sum(weights,ndim);
    stdev = sqrt(A ./ B);
else
    sq = (dat-me2).^2;
    sq(ib)=0; %set these to zero so they aren't included in the sum
    stdev = sqrt(  sum( sq, ndim ) ./ (nnums-1)   );
end



me(nnums==0)=NaN; %set the answer to NaN when have no non-NaN data for an entry
stdev(nnums==0)=NaN; %set the answer to NaN when have no non-NaN data for an entry
stdev(nnums-1==0)=0; %points where only had one value - can't divide by N-1 - set std to 0

if nargin>2
    if isstr(op)
        if strcmp(op,'sum')==1
            me=me.*nnums; %return the sum rather than the mean
        end
    end
end

 if isqueeze==1
    me=squeeze(me);  
    nnums=squeeze(nnums);
    stdev=squeeze(stdev);
 end
    








