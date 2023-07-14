function [nc_out]=ndhistc_run(varargin)
% Make sure have square brackets around the input data
% input format is [nc]=ndhistc_run([data01,data02....] ,bins01,bins02....)
% Wrapper to run ndhistc - cuts out data that is outside of the bins and NaN data
% as this can mess up ndhistc


[status,host]=system('echo $HOSTNAME');

ndims=nargin-1;

%max_chunk=2e7/ndims;
%nchunk = floor(size(varargin{1},1)/max_chunk);



icut=[];

for idims=1:ndims

    bins = varargin{1+idims}(:);
    %dat = varargin{1}(:,idims);

    ilow = find(varargin{1}(:,idims)<bins(1));
    ihigh = find(varargin{1}(:,idims)>bins(end));
    inan=find(isnan(varargin{1}(:,idims))==1);
    
    if length(ihigh)>0
        icut = unique([icut ihigh(:)']);
    end
    if length(ilow)>0
        icut = unique([icut ilow(:)']);
    end
    if length(inan)>0
        icut = unique([icut inan(:)']);
    end
end

%for idims=1:ndims
    varargin{1}(icut,:) = '';
%    varargin{(idims-1)*2+1}(length(dat)+1:end) = '';
%end


%for i=1:nchunk
    %nc=ndHistc(varargin{1:ndims},varargin{ndims+1:nargin});
%end

%Split the data into chunks and add the histograms at the end - trying to
%avoid the memory allocation issues, but it doesn't seem to be related to
%how big the data is...
ndat = size(varargin{1},1);
if ndims>2
    chunk=1;
    max_chunk=1e5;   
    nchunks = ceil(ndat / max_chunk)    
else
    chunk=0;
    nchunks=1;
end

istart=1;
for ichunk=1:nchunks
    ichunk
    if chunk==1
        iend = min(istart+max_chunk-1,ndat);
        inds = [istart:iend];
        istart=istart+max_chunk;
    else
        inds = [1:ndat];
    end
    
    new_dat_in{1} = varargin{1}(inds,:);
    for iarg=2:nargin
        new_dat_in{iarg}=varargin{iarg};
    end
    
    switch host(1:end-1)
        case 'challenger'
            nc=ndHistc_chall(varargin{1:ndims},varargin{ndims+1:nargin});
        otherwise
            %nc=ndhist(varargin{1:ndims},varargin{ndims+1:nargin}); %trying the matlab version due to issues with C version on Olympus
            %Takes longer, but is perhaps not prohibitive
            %nc=ndHistc(varargin{1:ndims},varargin{ndims+1:nargin});
            %nc=ndHistc_chall(varargin{1:ndims},varargin{ndims+1:nargin}); %sometimes this works on Olypmus
            %The non-C version has started being very slow - too many Zbins for
            %3D version?
            
            %chunked version
            % C-versions :-
            %nc=ndHistc_chall(new_dat_in{1:ndims},new_dat_in{ndims+1:nargin}); %sometimes this works on Olypmus
            %nc=ndHistc(new_dat_in{1:ndims},new_dat_in{ndims+1:nargin}); %sometimes this works on Olypmus
            % Matlab version :-
            nc=ndhist(new_dat_in{1:ndims},new_dat_in{ndims+1:nargin});
    end
    
    if ndims==1
        nc=nc(:,1); %If only have one dimension then it will still return an NxN matrix - just take the first column
    end
    
    if ichunk==1
        nc_out = zeros(size(nc));
    end
    
    nc_out = nc_out + nc;
    
end
