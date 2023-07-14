function X=Read1D2_v2(SIZE1,SIZE2,fid)
% Reads in a record from davie format
X=fread(fid,[SIZE1 SIZE2],'float=>double');

% Now to Skip the remainder of the record
REC_SIZE=0;
while(SIZE1.*SIZE2>REC_SIZE)
    REC_SIZE=REC_SIZE+128;
end
SKIP_SIZE=REC_SIZE-SIZE1.*SIZE2;
%fseek(fid,4.*SKIP_SIZE,'cof');
