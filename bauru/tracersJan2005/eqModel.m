clear dzz q dq qq dq_diff

qbase=5; %ppmv value that the values in dq_tot are based on


idir=1;

itdehyd=findheight(GridDan(1).t+3,23.75); %time index to take dq_tot values from

itdehyd=49;

z=GridDan(idir).Z+620;
     
[z0 zend]=findheight(z,15.8e3,17e3);

N=length(GridDan(idir).Y1);

rho=GridDan(idir).RHON(z0:zend);
zzeq=z(z0:zend);
dzz=(zzeq(2:end)-zzeq(1:end-1))/1000;
dzz(end+1)=dzz(end);

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
q(1)=5;
q(2:length(zzeq))=q(1); %starting values in ppmv - values given to the profile turn out to be irrelevant
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
     
            
    for iz=2:length(zzeq)
        
        rhoav=(rho(iz-1)+rho(iz))/2;
        
        dqdt_as(iz,t) = w*(rho(iz-1)*q(iz-1,t-1) - rho(iz)*q(iz,t-1) ) / rhoav / dzz(iz) ;
        
        dqdt_as(iz,t) = dqdt_as(iz,t) - w* (rho(iz-1) - rho(iz)) * q(iz)/dzz(iz);
        
        if igrad==1 %for gradual removal of vapour - effect spread evenly over time between events (1/N) dq/dt= change/dt = N*change
           % dqdt_ov(iz,t) =  Nev* prefact(iz)/rhoav/dzz(iz);  %spread based on the sum of dq (dq3) over all times. Done seperately for positive 
               %dqdt_ov(iz,t) = - Nev*dq_tot(idir).d(iz+z0-1,itdehyd,2);
            dqdt_ov(iz,t) = - Nev * ( dq_tot(idir).d(iz+z0-1,itdehyd,2) + nn(idir).n(iz+z0-1,2)/N *(q(iz)-qbase) ); %qbase=q deficit calc'd from e.g.=5ppmv 
                            %dq/dt=( Nev*dqtot/N + n/N(qenv-5) )    3rd index is 2 as refers to dqtot for <5 ppmv (see Allimp..)
                            %no need to divide dq_tot by length(Grid.Y1) since was done in Allimp..
        else
            dqdt_ov(iz,t)=0;
            
        end
        
        
        q(iz,t)=q(iz,t-1) + ( dqdt_as(iz,t) + dqdt_ov(iz,t) ) * dt;
        q(iz,t)=max( [0 q(iz,t)] ); %keep >=0
        
        if tovflag==1 & igrad==0
                %q(iz,t)=q(iz,t) + ( AA(iz) - BB(iz)*q(iz,t-1) ); %deplete q if time for overshoot
                q(iz,t)=max( [ q(iz,t) - dqz(iz)/sumdqz /rhoav/dzz(iz) 0 ] ); %keep above zero
        end
        
      
        
    end
    
%break    
end  %t=2:nt

%break

[qmin(n) imin]=min(q(2:end,end));
hmin(n)=zzeq(imin+1);
qq(:,n)=q(:,end);

end %n

'done'