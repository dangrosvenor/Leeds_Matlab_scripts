%3d version of dq dehyd in Allimp.. that reads info in gradually and saves e.g the total water to disk 
%in a seperate file

            '3d dq dehyd'
            f=1e6*28.97/18;
            imax=APARMS(1)+2;
            jmax=NJ+2;
            kmax=APARMS(3);
            
     for iq=1:6   
         if iq~=1
             fid2=fopen('totwater3d','rwb'); %open the file to store 3d tot water field
         end
             
         for k=1:kmax
             if iq~=1
                 %read in the running total of the kth slice
                 Xold=fread(fid2,[imax.*jmax],'float=>double');
			     Xold=reshape(Xold,imax,jmax); 
             else
                 Xold=zeros([imax jmax]);
             end
            
            %read in new value for iq
            X=fread(fid,[imax.*jmax],'float=>double');
			X=reshape(X,imax,jmax); 
            
            Xold=Xold+X;
            
            fseek(fid2,-imax*jmax,'cof'); %rewind the stored total for overwriting
            
            fwrite(fid2,Xold,'double');
        end
        fclose(fid2);
    end
            
    break
 
            
            j
            prcs=[0:5:100];
            tot_prctiles(j).t(1:size(totw,1),jj,1:length(prcs))=(prctile(totw',prcs))';
            vap_prctiles(j).t(1:size(vap,1),jj,1:length(prcs))=(prctile(vap',prcs))';
            
            minpps=[3.67 5 1 2 3 4]/f;
            %minpps=3.67/f;
          
            P=TwoD(idir).PP; %tot P
            thref=repmat(GridDan(idir).THREF,[1 length(GridDan(idir).Y1)]);
            T=TwoD(idir).TH1+thref; %tot potemp
            T=T./(1e5./P).^0.286; %tot temp                                                
            rho=P.*28.97e-3/8.3144./T;
            rhomoist=rho.*(1+TwoD(idir).Q(:,:,1))./(1+1.608*TwoD(idir).Q(:,:,1));  
  
          for ipps=1:length(minpps)
              
            for ikm=1:size(totw,1)

                                  
                    inon=find(totw(ikm,:)<minpps(ipps));
                    dq(j).d(ikm,jj,ipps)=sum(minpps(ipps)-totw(ikm,inon)) * f / size(totw,2);
                    nn(j).n(ikm,jj,ipps)=length(inon);

                    
                    inon2=find(vap(ikm,:)<minpps(ipps));
                    dq2(j).d(ikm,jj,ipps)=sum(minpps(ipps)-vap(ikm,inon2)) * f / size(vap,2);
                    nn2(j).n(ikm,jj,ipps)=length(inon2);
                    
                    if ipps==2
                        rho_prof(j).tot(ikm,jj)=mean(rhomoist(ikm,inon),2); %mean density in low tot water points
                        rho_prof(j).vap(ikm,jj)=mean(rhomoist(ikm,inon2),2); %same for low vapour points
                        rho_prof(j).mean(ikm,jj)=mean(rhomoist(ikm,:),2); %mean for all points
                    end
                    
              end    

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
            

