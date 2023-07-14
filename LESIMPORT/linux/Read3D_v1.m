function X=Read3D_v1(ISIZE,JSIZE,KSIZE,fid)

X=fread(fid,[ISIZE.*JSIZE.*KSIZE],'float=>double');
X=reshape(X,ISIZE,JSIZE,KSIZE);    
% Now to Skip the remainder of the record
REC_SIZE=0;
while(ISIZE.*JSIZE.*KSIZE>REC_SIZE)
    REC_SIZE=REC_SIZE+128;
end
SKIP_SIZE=REC_SIZE-ISIZE.*JSIZE.*KSIZE;
%fseek(fid,4.*SKIP_SIZE,'cof');
