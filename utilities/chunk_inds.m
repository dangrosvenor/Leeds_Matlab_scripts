function [start_inds,end_inds]=chunk_inds(N,maxN)

nchunks = ceil(N/maxN);
start_inds(1)=1;
for ichunk=1:nchunks
    end_inds(ichunk) = start_inds(ichunk)+maxN-1;
    if end_inds(ichunk)>=N
        end_inds(ichunk)=N;
        break
    else
        start_inds(ichunk+1) = end_inds(ichunk) + 1;
    end    
end

