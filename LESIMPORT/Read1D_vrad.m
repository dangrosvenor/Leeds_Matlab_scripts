function [X,part]=Read1D_vrad(SIZE,fid,linux,part)
% Reads in a record from davie format

if linux==0
    X=fread(fid,SIZE,'float=>double');
% else
% 
% 	%fseek(fid,SIZE,'cof');
% 	%X=ones([SIZE,1]);
% 	% Now to Skip the remainder of the record
% 	RECSIZE=128;
% 	REC_SIZE=0;
% 	irec=0;
% 	while(SIZE>=REC_SIZE)
%         temp(irec*RECSIZE+1:(irec+1)*RECSIZE)=fread(fid,RECSIZE,'float=>double');
%         REC_SIZE=REC_SIZE+RECSIZE;
%         irec=irec+1;
%         fseek(fid,4*RECSIZE,'cof');        
% 	
% 	end
% 	SKIP_SIZE=1024-SIZE;
% 	%fseek(fid,4.*SKIP_SIZE,'cof');
% 	
% 	X(1:SIZE)=temp(1:SIZE);
% 
% end

elseif linux==1
    % Reads in a record from davie format (Paul's routine)
	X=fread(fid,SIZE,'float=>double');
	
	% Now to Skip the remainder of the record
	REC_SIZE=0;
	while(SIZE>=REC_SIZE)
        REC_SIZE=REC_SIZE+128;
	end
	SKIP_SIZE=REC_SIZE-SIZE;
	%fseek(fid,4.*SKIP_SIZE,'cof');
    
elseif linux==2 %data is in blocks of 128 with a gap of 384 in between on Newton
    SKIP_SIZE=384;
    
    if SIZE>=part %if remainder of 128 block is to be read
        X(1:part)=fread(fid,part,'float=>double');
        n=part+1; %n is position of array to fill at next
        fseek(fid,4.*SKIP_SIZE,'cof');
    else
        X(1:SIZE)=fread(fid,SIZE,'float=>double');
        part=part-SIZE; %reduce the part of the 128 block that is left
        break; %and exit
    end

    rem=SIZE-part;
	while(rem>=128)

          X(n:n+127)=fread(fid,128,'float=>double');
          fseek(fid,4.*SKIP_SIZE,'cof');
          rem=rem-128;
          n=n+128;
        
    end
    
    if rem~=0
        X(n:n+rem-1)=fread(fid,rem,'float=>double'); %read in part of array remaining and leave file at this location
    end
    
	part=128-rem; %if rem==0 then part=128 and will read in first block in "if SIZE>=part" bit
    %will have skipped on ready for next read in while loop
	
    
end