clear dzz q dq qq dq_diff

qbase=5; %ppmv value that the values in dq_tot are based on


idir=1;

itdehyd=findheight(GridDan(1).t+3,23.75); %time index to take dq_tot values from

itdehyd=49;

N=length(GridDan(idir).Y1);


w=1.2; %ascent rate km/month

dt=1/60; %timesteps per month
tend=12; %in months
nt=round(tend/dt);


z=GridDan(idir).Z+620;     
[z0 zend]=findheight(z,15.8e3,17e3);

rho=GridDan(idir).RHON(z0:zend);
zz=z(z0:zend);
dzz=(zz(2:end)-zz(1:end-1))/1000;
dzz(end+1)=dzz(end);

nfi=5; %number of fine grid points moved up by slow ascent each time step
dzzf=w*dt/nfi; %fine grid spacing (km) calc'd so that are 10 fine levels covering air moved in each timestep
zzf=zz(1)/1000:dzzf:zz(end)/1000; %fine grid for vapour advection
rhof=interp1(zz/1000,rho,zzf); %interpolated density for fine grid
dqtotf=interp1(zz/1000,dq_tot(idir).d(z0:zend,itdehyd,2),zzf);
nntotf=interp1(zz/1000,nn(idir).n(z0:zend,itdehyd,2),zzf);

f=1e6*28.97/18;
T=tempLES(GridDan(idir)); %K
ei=SatVapPress(T(z0:zend),'goff','ice'); %Pa
P=GridDan(idir).PREFN(z0:zend); %Pa
sat=f*0.622*ei./(P-ei);




q(1)=4.5;
q(1)=5;
q(2:length(zzf))=q(1); %starting values in ppmv - values given to the profile turn out to be irrelevant
q=repmat(q,[nt 1])';

Nevs=[0:10:120]/10;

itend=36; %index of time to take cumulative sum over for overshoot mass deficit calc 

rho2=repmat(rho,[1 itend]);
dzz2=repmat(dzz,[1 itend]);


for n=1:length(Nevs)
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
     
    rhoav=(rhof(1:end-1)+rhof(2:end))/2;
    
     
    for iz=2:length(zzf)
        

        
        %dqdt_as(iz,t) = dqdt_as(iz,t) - w* (rho(iz-1) - rho(iz)) * q(iz)/dzz(iz);
        
 %       if igrad==1 %for gradual removal of vapour - effect spread evenly over time between events (1/N) dq/dt= change/dt = N*change
           % dqdt_ov(iz,t) =  Nev* prefact(iz)/rhoav/dzz(iz);  %spread based on the sum of dq (dq3) over all times. Done seperately for positive 
               %dqdt_ov(iz,t) = - Nev*dq_tot(idir).d(iz+z0-1,itdehyd,2);
%            dqdt_ov(iz,t) = - Nev * ( dq_tot(idir).d(iz+z0-1,itdehyd,2) + nn(idir).n(iz+z0-1,2)/N *(q(iz)-qbase) ); %qbase=q deficit calc'd from e.g.=5ppmv 
                            %dq/dt=( Nev*dqtot/N + n/N(qenv-5) )    3rd index is 2 as refers to dqtot for <5 ppmv (see Allimp..)
                            %no need to divide dq_tot by length(Grid.Y1) since was done in Allimp..
                            
                dqdt_ov = - Nev * ( dqtotf(iz) + nntotf(iz)/N *(q(iz,t)-qbase) ); %qbase=q deficit calc'd from e.g.=5ppmv 
                
                %       else
%            dqdt_ov(iz,t)=0;
            
%        end
        
        
        q(iz,t)=q(iz,t-1) + ( dqdt_ov ) * dt;
        q(iz,t)=max( [0 q(iz,t)] ); %keep >=0
        
%        if tovflag==1 & igrad==0
                %q(iz,t)=q(iz,t) + ( AA(iz) - BB(iz)*q(iz,t-1) ); %deplete q if time for overshoot
%                q(iz,t)=max( [ q(iz,t) - dqz(iz)/sumdqz /rhoav/dzz(iz) 0 ] ); %keep above zero
%        end
        
      
        
    end
    
    q(nfi+1:end,t)=q(1:end-nfi,t); %nfi is number of points shifted up each timestep
    q(1:nfi,t)=qbase;         %bring in qbase air from below 
    
%break    
end  %t=2:nt

%break

[qmin(n) imin]=min(q(2:end,end));
hmin(n)=zzf(imin+1);
qq(:,n)=q(:,end);

end %n

'done'

zzeq=zzf*1000; %for plotting in plotTime..