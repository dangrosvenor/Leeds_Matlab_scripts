function [X,part]=Read2D_v2(ISIZE,JSIZE,fid,linux,part)

if linux ~=2
    
X=fread(fid,[ISIZE,JSIZE],'float=>double')';

if linux==1
    rec=128;
else
    rec=512;
end

%fseek(fid,ISIZE*JSIZE,'cof');
%X=ones([ISIZE,JSIZE]);
% Now to Skip the remainder of the record
REC_SIZE=0;
while(ISIZE.*JSIZE>REC_SIZE)
    REC_SIZE=REC_SIZE+rec;
end
SKIP_SIZE=REC_SIZE-ISIZE.*JSIZE;
%fseek(fid,4.*SKIP_SIZE,'cof');

%method for multi-d arrays is to read in as a 1-d array and reshape at end
elseif linux==2 %data is in blocks of 128 with a gap of 384 in between on Newton
    SKIP_SIZE=384;
    SIZE=ISIZE*JSIZE;
    
    if SIZE>=part %if remainder of 128 block is to be read
        x(1:part)=fread(fid,part,'float=>double');
        n=part+1; %n is position of array to fill at next
        fseek(fid,4.*SKIP_SIZE,'cof');
    else
        x(1:SIZE)=fread(fid,SIZE,'float=>double');
        part=part-SIZE; %reduce the part of the 128 block that is left
        break; %and exit
    end

    rem=SIZE-part;
	while(rem>=128)

          x(n:n+127)=fread(fid,128,'float=>double');
          fseek(fid,4.*SKIP_SIZE,'cof');
          rem=rem-128;
          n=n+128;
        
    end
    
    if rem~=0
        x(n:n+rem-1)=fread(fid,rem,'float=>double'); %read in part of array remaining and leave file at this location
    end
    
	part=128-rem; %if rem==0 then part=128 and will read in first block in "if SIZE>=part" bit
    %will have skipped on ready for next read in while loop
    
    X=reshape(x,[ISIZE JSIZE]);
	
    
end
