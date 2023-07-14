function X=Read1D_vrad(SIZE,fid)
% Reads in a record from davie format
X=fread(fid,SIZE,'float=>double');

% Now to Skip the remainder of the record
REC_SIZE=0;
while(SIZE>=REC_SIZE)
    REC_SIZE=REC_SIZE+1024;
end
SKIP_SIZE=REC_SIZE-SIZE;
fseek(fid,4.*SKIP_SIZE,'cof');
