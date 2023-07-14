function X=Read1D2_vrad(SIZE1,SIZE2,fid,linux)
% Reads in a record from davie format
X=fread(fid,[SIZE1 SIZE2],'float=>double');

if linux==1
    rec=128;
else
    rec=512;
end

% Now to Skip the remainder of the record
REC_SIZE=0;
while(SIZE1.*SIZE2>REC_SIZE)
    REC_SIZE=REC_SIZE+rec;
end
SKIP_SIZE=REC_SIZE-SIZE1.*SIZE2;
%fseek(fid,4.*SKIP_SIZE,'cof');


