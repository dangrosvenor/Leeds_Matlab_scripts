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
    
    N=floor( (SIZE-part)/128 ); %no. of 128 sized blocks left 
    
    if SIZE>=part %if remainder of 128 block is to be read
        X(1:part)=fread(fid,part,'float=>double');
        n=part+1; %n is position of array to fill at next
        fseek(fid,4.*SKIP_SIZE,'cof');
    else
        X(1:SIZE)=fread(fid,SIZE,'float=>double');
        part=part-SIZE; %reduce the part of the 128 block that is left
        X=X';
        break; %and exit
    end
    
    x=fread(fid,512*N,'float=>double'); %read in all N blocks with gaps in between
    
%     for j=1:N
%         ind((j-1)*384+1:j*384)=[(j-1)*512+129:j*512]; %declare indices where the gaps are
%     end
    
    b=repmat([129:512],[N 1]); %make matrix of 129:512 in each row for N rows
    d=repmat([1:N]*512,[384 1]); %make matrix of 512s, 1024s, etc for each row of N
    d=d-512; %take away one 512 as want first row to be 0s, 2nd 512s, etc.
    d=d';
    e=b+d; %add (i-1)*512 to each row, i, of N
    ind=reshape(e,[1 N*384]); %unroll 2d matrix into vector
    %ind now = matrix of the indices of all the gaps
    
    if exist('ind')
        x(ind)=[]; %remove the gaps from x
    end
    
    X(n:n+N*128-1)=x(1:length(x));
    
    rem=SIZE-part-N*128;
    
    n=n+N*128;
    
    
%     rem=SIZE-part;
% 	while(rem>=128)
% 
%           X(n:n+127)=fread(fid,128,'float=>double');
%           fseek(fid,4.*SKIP_SIZE,'cof');
%           rem=rem-128;
%           n=n+128;
%         
%     end
%     
    if rem~=0
        X(n:n+rem-1)=fread(fid,rem,'float=>double'); %read in part of array remaining and leave file at this location
    end
    


	part=128-rem; %if rem==0 then part=128 and will read in first block in "if SIZE>=part" bit
    %will have skipped on ready for next read in while loop
    
    X=reshape(X,[ISIZE JSIZE]);
    
    X=X';
	
    
end
