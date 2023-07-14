function X=Read3D(imax,jmax,kmax,fid)


for k=1:kmax
   
    %fseek(fid,8,'cof');
    for j=1:jmax
        if(j==2)
            %fseek(fid,8,'cof');
        end
        X(:,j,k)=fread(fid,imax,'float=>double');
        if(j==jmax-1)
           %fseek(fid,8,'cof');
        end
        
    end
    for j=1:jmax
        if (k<=3)
       fseek(fid,8,'cof');
       
            fseek(fid,8,'cof');
       end
   end
end

