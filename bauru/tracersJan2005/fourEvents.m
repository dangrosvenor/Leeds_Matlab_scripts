clear dzz q


idir=2;
iev=[11 14 20 24]; %time indices of 4 events

dmcum=1.8424/1000; %sum of positive increases only (kg)

Am=3.67 - 502./nn(2).n(:,iev) .* dq_tot(idir).d(:,iev,1);
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
q(2:length(zz))=q(1); %starting values in ppmv - values given to the profile turn out to be irrelevant
q=repmat(q,[nt 1])';

Nevs=[10:10:120];


for n=1:length(Nevs)

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
            
    for iz=2:length(zz)
        rhoav=(rho(iz-1)+rho(iz))/2;
        dqdt_as(iz,t) = w*(rho(iz-1)*q(iz-1,t-1) - rho(iz)*q(iz,t-1) ) / rhoav / dzz(iz) ;
        
        
        
        if igrad==1
            %dqdt_ov(iz,t) = Nev* ( AA(iz) - BB(iz)*q(iz,t-1) );    
            dqdt_ov(iz,t) = - Nev* dmcum/(length(zz)-1) /rhoav/dzz(iz);
        else
            dqdt_ov(iz,t)=0;
            
        end
        
        
        q(iz,t)=q(iz,t-1) + ( dqdt_as(iz,t) + dqdt_ov(iz,t) ) * dt;
        
        if tovflag==1 & igrad==0
                q(iz,t)=q(iz,t) + ( AA(iz) - BB(iz)*q(iz,t-1) ); %deplete q if time for overshoot
        end
        
    end
    
end

'done';
qmin(n)=min(q(2:end,end));

end

'done'