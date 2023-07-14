function [X,part]=Read2D_v2(ISIZE,JSIZE,fid,linux,part)

if linux~=1
X=fread(fid,[ISIZE,JSIZE],'float=>double')';
% Now to Skip the remainder of the record
%REC_SIZE=0;
%while(ISIZE.*JSIZE>REC_SIZE)
%    REC_SIZE=REC_SIZE+512;
%end
%SKIP_SIZE=REC_SIZE-ISIZE.*JSIZE;
%fseek(fid,4.*SKIP_SIZE,'cof');

else
    
X=fread(fid,[ISIZE,JSIZE],'float=>double')';
% Now to Skip the remainder of the record
REC_SIZE=0;
while(ISIZE.*JSIZE>REC_SIZE)
    REC_SIZE=REC_SIZE+128;
end
SKIP_SIZE=REC_SIZE-ISIZE.*JSIZE;
%fseek(fid,4.*SKIP_SIZE,'cof');

end