%Run Diag_3d_newt with idontclose=1 to stop the import just before the Q-fields are due to be read in
clear horiz

fprintf(1,'\n Loading file...');
Diag_3d_newt;
%GridDan(1)=Grid;

            fprintf(1,'\n3d slice getter running.....');
            f=1e6*28.97/18;
            imax=APARMS(1)+2;
            jmax=NJ+2;
            kmax=APARMS(3);
 


izpos=99;

%ixpos=2;
%iypos=2;
 
        
     for iq=1:9   
  
         for k=1:kmax             
            
            %read in new value for iq
            X=fread(fid,[imax.*jmax],'float=>double');           
            X=reshape(X,[imax jmax]);
            X=X(2:end-1,:);
           
            if k==izpos
               horiz(:,:,iq)=X;
           end
		%	TwoD.Q(k,:,iq)=X(:,ixpos)';                                  
        end
        
    end
                    
                    
		
            
    fprintf(1,'\ndone');
%    break
    
