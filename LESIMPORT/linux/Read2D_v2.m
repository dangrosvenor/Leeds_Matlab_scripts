function X=Read2D_v2(ISIZE,JSIZE,fid)

X=fread(fid,[ISIZE,JSIZE],'float=>double')';
% Now to Skip the remainder of the record
REC_SIZE=0;
while(ISIZE.*JSIZE>REC_SIZE)
    REC_SIZE=REC_SIZE+128;
end
SKIP_SIZE=REC_SIZE-ISIZE.*JSIZE;
%fseek(fid,4.*SKIP_SIZE,'cof');
