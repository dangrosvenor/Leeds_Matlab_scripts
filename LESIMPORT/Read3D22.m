function X=Read3D(imax,jmax,kmax,fid)

% for i=1:2
% for j=1:jmax*imax
%    fseek(fid,8,'cof');
% end
% end

for k=1:kmax
    for j=1:jmax
        
        if (k<=1)
             %fseek(fid,16,'cof');
%             fseek(fid,8,'cof');   
       end
   end
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
    
end

