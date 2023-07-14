function [dat_out,time_out] = sort_multidim_array_time(dat,time_glm,chunk_dim)
%function [dat_out,time_out] = sort_multidim_array_time(dat,time_glm,chunk_dim)
%Sorts the dat array according to time_glm, assuming that the time
%dimension for dat is the 1st dimension.
%Allows the data to be processed in batches along chunk_dim (not the time
%dimension) in case it runs into memory issues.

sdat = size(dat);

[time_out,I]=sort(time_glm);

if ~exist('chunk_dim')
    chunk_dim=0;
end

if chunk_dim==0
    chunk_str = ':,';
    Nchunk = 1;
else
    chunk_str = 'ichunk,';
    Nchunk = sdat(chunk_dim);
end

for ichunk=1:Nchunk
    %replace this with version where can choose chunk_dim
    %[dat_chunk] = sort_multidim_chunk(dat(:,ichunk,:,:,:),I);
    eval_str=['[dat_chunk] = sort_multidim_chunk(dat(:,'];
    for i=2:length(sdat)
        if i==chunk_dim
            eval_str = [eval_str chunk_str];
        else
            eval_str = [eval_str ':,'];
        end
    end
    eval_str(end)=[];
    eval_str = [eval_str '),I);'];
    eval(eval_str);
    
    if ichunk==1
        dat_out = dat_chunk;
    else
        dat_out = cat(chunk_dim,dat_out,dat_chunk);
    end
end



function [dat_out] = sort_multidim_chunk(dat,I)

sdat = size(dat);
%sdat = 1980         144         192  [time,lat,lon]
rep_inds=1;
for idim=2:length(sdat)
    rep_inds=[rep_inds sdat(idim)];
end

N = prod(sdat);


%Get subscripts for the whole 3d array
%Replace the following with an N-dim version
%[i1,i2,i3] = ind2sub(sdat,[1:prod(sdat)]);
eval_str='[';
for i=1:length(sdat)
    eval_str = [eval_str 'i' num2str(i) ','];
end

eval_str(end)=[]; %remove the extra comma
eval_str = [eval_str ']=ind2sub(sdat,[1:N]);'];
eval(eval_str);





%New indices for just the time dim using I values for
%sorting - replicate to 3D array
i1_new = repmat(I,rep_inds);

%Get indices using new values
%inew = sub2ind(sdat,i1_new(:)',i2,i3);
%convert to multi-dim version :-
eval_str='inew = sub2ind(sdat,i1_new(:)''';
for i=2:length(sdat)
    eval_str = [eval_str ',i' num2str(i)];
end
eval_str = [eval_str ');'];
eval(eval_str);

%Use to sort dat and re-shape
dat_out = dat(inew);
dat_out = reshape(dat_out,sdat);





%
%
% N = prod(sdat);
% maxN = 1e8;
% nchunks = ceil(N/maxN);
% istart=1;
% for ichunk=1:nchunks
%     iend = istart+maxN-1;
%     if iend>=N
%         iend=N;
%         last=1;
%     end
%     sdat_chunk = sdat;
%     sdat_chunk(1) = length([istart:iend]);
%
%     %Get subscripts for the whole 3d array
%     %Replace the following with an N-dim version
%     %[i1,i2,i3] = ind2sub(sdat,[1:prod(sdat)]);
%     eval_str='[';
%     for i=1:length(sdat)
%         eval_str = [eval_str 'i' num2str(i) ','];
%     end
%
%     eval_str(end)=[]; %remove the extra comma
%     eval_str = [eval_str ']=ind2sub(sdat,[istart:iend]);'];
%     eval(eval_str);
%
%
%
%
%
%     %New indices for just the time dim using I values for
%     %sorting - replicate to 3D array
%     i1_new = repmat(I,rep_inds);
%
%     %Get indices using new values
%     %inew = sub2ind(sdat,i1_new(:)',i2,i3);
%     %convert to multi-dim version :-
%     eval_str='inew = sub2ind(sdat,i1_new(:)''';
%     for i=2:length(sdat)
%         eval_str = [eval_str ',i' num2str(i)];
%     end
%     eval_str = [eval_str ');'];
%     eval(eval_str);
%
%     %Use to sort dat and re-shape
%     dat_out = dat(inew);
%     dat_out = reshape(dat_out,sdat);
%
%     if last==1
%         break
%     else
%         istart = iend + 1;
%     end
%
% end
%
%
%
%
%


        