clear dzz q dq qq dq_diff

qbase=5; %ppmv value that the values in dq_tot are based on


idir=1;
iev=[11 14 20 24]; %time indices of 4 events

dmcum=1.8424/1000; %sum of positive increases only (kg)

Am=3.67 - 502./nn(idir).n(:,iev) .* dq_tot(idir).d(:,iev,1);
Am(isnan(Am))=0;
A=(sum(nn(idir).n(:,iev).*Am,2))/502;
B=(sum(nn(idir).n(:,iev),2))/502;


z=GridDan(idir).Z+620;
     
[z0 zend]=findheight(z,15.8e3,17e3);


rho=GridDan(idir).RHON(z0:zend);
zz=z(z0:zend);
dzz=(zz(2:end)-zz(1:end-1))/1000;
dzz(end+1)=dzz(end);
AA=A(z0:zend);
BB=B(z0:zend);


f=1e6*28.97/18;
T=tempLES(GridDan(idir)); %K
ei=SatVapPress(T(z0:zend),'goff','ice'); %Pa
P=GridDan(idir).PREFN(z0:zend); %Pa
sat=f*0.622*ei./(P-ei);



w=1.2; %ascent rate km/month

dt=1/600; %timesteps per month
tend=6; %in months
nt=round(tend/dt);

q(1)=4.5;
q(2:length(zz))=4.5; %starting values in ppmv - values given to the profile turn out to be irrelevant
q=repmat(q,[nt 1])';

Nevs=[10:10:120];

itend=36; %index of time to take cumulative sum over for overshoot mass deficit calc 

rho2=repmat(rho,[1 itend]);
dzz2=repmat(dzz,[1 itend]);



for iz=1:length(zz)
	dqz(iz)=sum(dq_tot(idir).d(iz+z0-1,1:itend,1),2);
    dqz2(iz)=sum(dq_tot(idir).d(iz+z0-1,1:itend,1)*rho(iz)*dzz(iz),2) ;
    dmass(iz)=rho(iz)*dzz(iz) ;
end
sumdqz=sum(dqz(2:end));
sumdqz2=sum(dqz2(2:end));

summass=sum(dmass(2:end));

%dqz2(:)=sumdqz2/(length(dqz2)-1);
%dqz2(3)=sumdqz2;
dqz2=dmass;
sumdqz2=summass;


for iz=1:length(zz)
    [maxdq(iz) imax]=max(dq_tot(idir).d(iz+z0-1,1:itend,1),[],2);
    nnmax(iz)=nn(idir).n(iz+z0-1,imax);
end

summax=sum(maxdq(2:end));


 for iz=1:length(zz)
         dq_st(iz,1:itend)=dq_tot(idir).d(iz+z0-1,1:itend,1); % - nn(idir).n(iz+z0-1,1:itend)/502 .* (3.67 - q(iz,t-1) );
         %dq(iz,1:itend)=dq_tot(idir).d(iz+z0-1,1:itend,1) - nn(idir).n(iz+z0-1,1:itend)/502 .* (3.67 - 3.8 );  
         
 end
     
     mtim=sum(dq_st.*rho2.*dzz2,1);
     
     dw=mtim(2:end)-mtim(1:end-1);
	 idw_neg=find(dw<0);
     idw_pos=find(dw>=0);





for n=1:12 %length(Nevs)
fprintf(1,'%d ',n);

Nev=Nevs(n); %number of events per month



igrad=1;
tov=1/Nev;

tt=0;
for t=2:nt
    tt=tt+dt;
    tovflag=0;
    
     if tt>=tov
        tovflag=1; %flag to say it's time to deplete vapour
        tt=0;
     end
     
  
     
     for iz=1:length(zz)
         dq(iz,1:itend)=dq_tot(idir).d(iz+z0-1,1:itend,1) - nn(idir).n(iz+z0-1,1:itend)/502 .* (qbase - q(iz,t-1) );
         %dq(iz,1:itend)=dq_tot(idir).d(iz+z0-1,1:itend,1) - nn(idir).n(iz+z0-1,1:itend)/502 .* (3.67 - 3.8 );  
          dq3(iz)=sum(dq(iz,idw_neg).*rho(iz).*dzz(iz),2); %measure of how long and how large dry points were over all times
                %for sum over all times to see if pos or neg - for scaling of dmcum - neg points based on how much negative influence overall
                %not just for contribution to dmcum - no points and magnitude

         
     end
     
     mtim=sum(dq.*rho2.*dzz2,1);
     
     dw=mtim(2:end)-mtim(1:end-1);
	 idw=find(dw<0);
	 %dw(idw)=0;
	 sw=cumsum(dw);
     
     %dmcum=sw(end);
     
     dw_new=dw(idw_neg); %dw for points where dissipation was deemed to be occuring from inital state (sum over time)
     dmcum=sum(dw_new); %using positve or negative dm changes from original values based on 3.67 ppmv
     
     idw_neg3=find(dw_new<0);
     idw_pos3=find(dw_new>0);
     
     dmcum_neg=sum(dw_new(idw_neg3));
     dmcum_pos=sum(dw_new(idw_pos3));
     
     ipos3=find(dq3>0);
     ineg3=find(dq3<0);
     
     dq3tot_pos=sum(dq3(ipos3)); %positive sign
     dq3tot_neg=sum(dq3(ineg3)); %negative sign

          
     
     
     prefact(ipos3)=dmcum_neg*dq3(ipos3)/dq3tot_pos;
     prefact(ineg3)=dmcum_pos*dq3(ineg3)/dq3tot_neg;
     
     if (length(ineg3)>=1)
         '';
     end
            
    for iz=2:length(zz)
        
        %dmcum pos for drying so prefact=pos for drying and neg for moistening
        %pos dq3 for moistening, neg for drying        
%         if length(find(ipos3==iz))>=1 
%             prefact(iz)=dmcum*dq3(iz)/dq3tot_pos; %drying pos prefact but uses neg dmcum_neg since neg of this means mixing with dry points
%         elseif length(find(ineg3==iz))>=1 
%             prefact(iz)=dmcum*dq3(iz)/dq3tot_neg; %moistening neg prefact
%         else 
%             prefact(iz)=0;
%         end
        
        
        rhoav=(rho(iz-1)+rho(iz))/2;
        
        dqdt_as(iz,t) = w*(rho(iz-1)*q(iz-1,t-1) - rho(iz)*q(iz,t-1) ) / rhoav / dzz(iz) ;
        
        dq_diff(iz,1:itend-1)=dq(iz,2:itend)-dq(iz,1:itend-1);
        %dq_diff=dq(iz,2:itend)-dq(iz,1:itend-1);
        
        if igrad==1
            %dqdt_ov(iz,t) = Nev* ( AA(iz) - BB(iz)*q(iz,t-1) );    
            %dqdt_ov(iz,t) = - Nev* dmcum/(length(zz)-1) /rhoav/dzz(iz);
            %dqdt_ov(iz,t) = - Nev* dmcum*dqz(iz)/sumdqz /rhoav/dzz(iz);
            %dqdt_ov(iz,t) = - Nev* (3*maxdq(iz) - nnmax(iz)/502 .* (3.67 - q(iz,t-1) )  );  
       %dqdt_ov(iz,t) =  Nev* dmcum*dqz2(iz)/sumdqz2 /rhoav/dzz(iz);  %spread based on mass
            %dqdt_ov(iz,t) = - Nev* dmcum*maxdq(iz)/summax /rhoav/dzz(iz);
            %dqdt_ov(iz,t) = 0;
            %dqdt_ov(iz,t) = - Nev * sum(dq_diff(idw_pos)); %sum over those times where original dm was increasing

            %dqdt_ov(iz,t) = - Nev * sum(dq_diff(iz,idw_pos),2); %sum over those times where original dm was increasing
            
            dqdt_ov(iz,t) =  Nev* prefact(iz)/rhoav/dzz(iz);  %spread based on the sum of dq (dq3) over all times. Done seperately for positive 
            %and negative heights - e.g. dmcum_neg*dq3(z)/dq3tot_pos (if positive) -graph name dq_posnegspread

        else
            dqdt_ov(iz,t)=0;
            
        end
        
        
        q(iz,t)=q(iz,t-1) + ( dqdt_as(iz,t) + dqdt_ov(iz,t) ) * dt;
        
        if tovflag==1 & igrad==0
                %q(iz,t)=q(iz,t) + ( AA(iz) - BB(iz)*q(iz,t-1) ); %deplete q if time for overshoot
                q(iz,t)=q(iz,t) - dqz(iz)/sumdqz /rhoav/dzz(iz);
        end
        
      
        
    end
    
%break    
end  %t=2:nt

%break

[qmin(n) imin]=min(q(2:end,end));
hmin(n)=zz(imin+1);
qq(:,n)=q(:,end);

end %n

'done