function [X, part, iread]=Read3D(imax,jmax,kmax,fid,linux,part,iread,nreads)

i=find(nreads==1);
if iread>max(i); return; end;

if nreads(iread)==1

	if linux~=1
        
		for k=1:kmax
            X(:,:,k)=fread(fid,[imax,jmax],'float=>double')';
		end
        
	else 

		X=fread(fid,[imax.*jmax.*kmax],'float=>double');
		X=reshape(X,imax,jmax,kmax);  
	%    X=reshape(X,JSIZE,ISIZE,KSIZE); 
	%  for k=1:250   
	%     temp=fread(fid,[ISIZE.*JSIZE],'float=>double');
	%     X(1:514,1:34,k)=reshape(temp,ISIZE,JSIZE);
	% end
		% Now to Skip the remainder of the record
		REC_SIZE=0;
		while(imax.*jmax.*kmax>REC_SIZE)
            REC_SIZE=REC_SIZE+128;
		end
		SKIP_SIZE=REC_SIZE-imax.*jmax.*kmax;
		%fseek(fid,4.*SKIP_SIZE,'cof');
	
	end
    
else
    for k=1:kmax
        fread(fid,[imax.*jmax],'float=>double'); %only read in cross sections to use less memory
    end
    %fseek(fid,4*(imax)*(jmax)*kmax,'cof'); %fseek doesn't seem to work!
    %clear X
end

iread=iread+1;