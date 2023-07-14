function X=Read1D_v3(SIZE,fid)
% Reads in a record from davie format
X=fread(fid,SIZE,'float=>double');

% Now to Skip the remainder of the record
