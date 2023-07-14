%3d version of dq dehyd in Allimp.. that reads info in gradually and saves e.g the total water to disk 
%in a seperate file
%set idontclose==1 in DIAG...

            '3d dq dehyd'
            f=1e6*28.97/18;
            imax=APARMS(1)+2;
            jmax=NJ+2;
            kmax=APARMS(3);
 
            
diag='th';
diag='vap';
diag='totw';

switch diag
case 'vap'
    D=14e3;
    [iy iy2]=findheight( Grid.Y1,Grid.Y1(1)+D/2,Grid.Y1(end)-D/2 );
    [ix ix2]=findheight( Grid.X1,Grid.X1(1)+D/2,Grid.X1(end)-D/2 );
    [iz2]=findheight( Grid.Z , 2.5e3 );
    
    vap(:,:,:,jj)=ThreeD.Q([iy2-1:imax 2:iy+2],[ix2-2:jmax 1:ix+2] , 1:iz2);
    
case 'th'
    D=14e3;
    [iy iy2]=findheight( Grid.Y1,Grid.Y1(1)+D/2,Grid.Y1(end)-D/2 );
    [ix ix2]=findheight( Grid.X1,Grid.X1(1)+D/2,Grid.X1(end)-D/2 );
    [iz2]=findheight( Grid.Z , 2.5e3 );
    
    theta(:,:,:,jj)=ThreeD.TH1([iy2-1:imax 2:iy+2],[ix2-2:jmax 1:ix+2] , 1:iz2);
    
case 'totw'
    
maketotw=0; %set this to zero and justvap to one for just the icediags_5thSept part
justvap=1; %flag to say to just do the vapour stats

prcs=[0:5:100];
minpps=[3.67 5 1 2 3 4]/f;
            
 nks=3;    
 kchunk=floor(kmax/nks);  
 ks=ones([1 nks-1])*kchunk;
 ks=[ks kmax-(nks-1)*kchunk]; %array of sizes of chunks to read in
%ks=[kchunk:kchunk:kchunk*(nks-1) kmax];
 
if justvap==1
    iqend=1;
else
    iqend=6;
end

if maketotw==1            
     for iq=1:iqend   

             fid2=fopen('totwater3d_w','wb'); %open the file to store 3d tot water field
             fid3=fopen('totwater3d_r','rb'); %open the file to read in old 3d tot water field
             
             if iq==1
                 fclose(fid3);  %close the file (opended and closed in order to make inactive) 
             end

  
         for k=1:kmax
             if iq~=1
                 %read in the running total of the kth slice
                 Xold=fread(fid3,[imax.*jmax],'double');
			   %  Xold=reshape(Xold,imax,jmax); 
             else
                 Xold=zeros([imax*jmax 1]);
                 
             end
            
            %read in new value for iq
            X=fread(fid,[imax.*jmax],'float=>double');

            Xold=Xold+X;
            
            X=reshape(X,[imax jmax]);
            X=X(2:end-1,:);
            X=reshape(X,[(imax-2)*jmax 1]);

            if iq==1  %if is vapour then do vapour diags here
                vap_prctiles(j).t(k,jj,1:length(prcs))=(prctile(X,prcs))';
                
                for ipps=1:length(minpps)
                    
                    inon2=find(X<minpps(ipps));
                    dq_vaps(j).d(k,jj,ipps)=sum(minpps(ipps)-X(inon2)) * f / length(Grid.Y1);
                    nn2(j).n(k,jj,ipps)=length(inon2);                                                              
                    
                end
                
                bins=[0:0.1:30];
                %dat=ThreeD.Q(2:end-1,:,ih1:ih2)*f;

                    vapdist(1).v(:,jj,k)=binner(f*X,bins); 
                    meanvap(1).m(:,jj,k)=mean(X);

   
                
                
            end
                    
                    
		%	X=reshape(X,imax,jmax); 
                    
            
            
            
            
            fwrite(fid2,Xold,'double');
        end
        
        fclose(fid2);
        if iq~=1
            fclose(fid3);
        end
        eval('!c:/cygwin/bin/mv totwater3d_w totwater3d_r'); %rename updated file for reading in next time

    end
    
    fclose(fid); %close diag file
            
    'done'
%    break
    
end
 
if justvap==1
    icediags_5thSept_2005_32;
    return  %exit from this routine and return to calling routine
end
        
%%% now do diags on saved data          
            
            fid=fopen('totwater3d_r','rb'); %open the file to store 3d tot water field
              
            for ikm=1:kmax
                %read in the kth slice
                 X=fread(fid,[imax.*jmax],'double');
			     X=reshape(X,imax,jmax); 
                 X=X(2:end-1,:);
                 X=reshape(X,[(imax-2)*jmax 1]);
                 
                 tot_prctiles(j).t(ikm,jj,1:length(prcs))=(prctile(X(2:end-1,:),prcs));
                
              for ipps=1:length(minpps)
                               
                    inon=find(X<minpps(ipps));
                    dq_tot(j).d(ikm,jj,ipps)=sum(minpps(ipps)-X(inon)) * f / length(Grid.Y1);
                    nn(j).n(ikm,jj,ipps)=length(inon);                                                  
                    
              end    
              
                    bins=[0:0.1:200];
                    totdist(1).v(:,jj,k)=binner(f*X,bins); 
                    meantot(1).m(:,jj,k)=mean(X);

            end
             
             icediags_5thSept_2005_32;

             
             %if jj>=8
%                 f=1e6*28.97/18;
%                 timei=mod(SER(end,1)/3600 + 19.67,24);
%                 timlab=num2str(round2(timei,2));
%                 noplot=1;
%                 plotTimeHeightVap3;
%                 output=1;
%                 exdirA=[direcDan(jc).dir 'results/vapMR+iceNC/'];
                %end
            

end