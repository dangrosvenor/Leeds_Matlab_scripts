%like eqfine but not storing q values at each time so saves memory and quicker
clear dzz q dq qq dq_diff qsave zw wref

qbase=5; %ppmv value that the values in dq_tot are based on (as used in Allimp.. - don't change unless using different deficit values
qsat=5; %MR of air entering from base of region (qsat for assumed cold point) %cold point mr of sounding =6.43

isatenv=1; %makes ice sat profile used equal to that of vapour profile from sounding (zero for ice sat at base, one for 5 ppmv)
isat=0; %flag to make the initial base input value the saturation MR of that height (one for ice sat at base)
iqenv=1; %flag to make it the environmental vapour MR  (zero for ice sat at base one for 5 ppmv)
ibelsat=1; %flag to keep below ice saturation (or below vapour profile)
ibelsatice=0; %flag to keep ice values used for moistening below ice sat
in2EQn1=1; %flag to make the number of points in the third dimension (n2) equal to the number in the 2-d sim (n)

imoisten=0;
Nmoist=1; %number of moistening events per month
idirCCN=2; %directory number for run to be used for moistening from ice

oldadvect=0; %flag to use old advection scheme where layers moved up (not really suitable for varying ascent rate with height)

wfactor=1; %factor to multiply ascent rates by for sensitivities

areafactor=0.5; %factor to multiply nntotf by to represent the dehydrated points covering less area in 3-d than the circular estimate

if ~exist('istore'); istore=0; end

H1=15.8e3;
H1=14e3;
H2=20e3; 
H1=14.8e3;


h0=20; %height above which to ignore dqtot and nn data
h0=17;

nclouds=1; %factor to account for no. of clouds like those represented in dqtot present at same time (more effect than if just increase event no. ...
            %due to it making so is less space for environmental air

N=length(GridDan(idir).Y1);
n2=N/5; %number of domain points affected in 3rd dimension



idir=1;

itdehyd=findheight(GridDan(1).t+3,23.75); %time index to take dq_tot values from
itdehyd=findheight(GridDan(1).t+3,23.84); %time of lowest mean tot water in 15.8 - 16 km region

itdehyd=62; %final profile from 1 km run
%itdehyd=36;



w=1.2; %ascent rate km/month

dt=1/360; %reciprocal of timesteps per month
%tend=18; %in months

%%%%%%%%%%%%%%%%time length of simulation
tend=15.5/30; %running for 8 days to see dehydration effect on horizontally moving parcel through cold pool region
            %assumes a region of about 3500km in size with air travelling at 5 m/s
nt=round(tend/dt);


clear diff
x2=n2*diff(GridDan(idir).Y1(1:2))/1000; %length in 3rd dimension (for plotTimeH....)

z=GridDan(idir).Z+620;   


[z0 zend]=findheight(z,H1,H2);

rho=GridDan(idir).RHON(z0:zend);
zz=z(z0:zend);
dzz=(zz(2:end)-zz(1:end-1))/1000;
dzz(end+1)=dzz(end);

nfi=5; %number of fine grid points moved up by slow ascent each time step

dzzf=w*dt/nfi; %fine grid spacing (km) calc'd so that are 10 fine levels covering air moved in each timestep
dzzf=0.1;
zzf=zz(1)/1000:dzzf:zz(end)/1000; %fine grid for vapour advection
rhof=interp1(zz/1000,rho,zzf); %interpolated density for fine grid
dqtotf=nclouds*interp1(zz/1000,dq_tot(idir).d(z0:zend,itdehyd,2),zzf);
nntotf=nclouds*interp1(zz/1000,nn(idir).n(z0:zend,itdehyd,2),zzf);

nntotf=nntotf*areafactor;

%set up reference w profile (from Jensen, 2001 paper - in balance with radiative heating. wref in cm/s
zw(1)=14;
wref(1)=-0.08;
zw(2)=16.6;
wref(2)=0.22;
zw(3)=16.8;
wref(3)=0.1; %0.1
zw(4)=25;
wref(4)=0.1;

% zw(4)=17;
% wref(4)=0.05;
% zw(5)=19;
% wref(5)=-0.01;

wref=wref*1e-5 * 3600*24*30 * wfactor; %convert from cm/s to km/month

ww=interp1(zw,wref,zzf); %interpolated velocity for fine grid

ww(ww<0)=0;


if in2EQn1==1
	n2=nntotf;
    n2str='n1';
else
    n2str=num2str(n2);    
    n2=repmat(n2,[1 size(nntotf,2)]);
end


ih0=findheight(zzf,h0);
dqtotf(ih0:end)=0;
nntotf(ih0:end)=0;

f=1e6*28.97/18;
T=tempLES(GridDan(idir)); %K
ei=SatVapPress(T(z0:zend),'goff','ice'); %Pa
P=GridDan(idir).PREFN(z0:zend); %Pa
sat=f*0.622*ei./(P-ei);
sat=interp1(zz/1000,sat,zzf); %interpolate onto fine grid

if isatenv==1
	sat=interp1(zz/1000,f*icediagsALL(idir).i(z0:zend,1,37)/npess2(idir),zzf); %if want to make ice sat MR = inital environmental MR of sounding
    sat(zzf>18.5)=5;
end

if imoisten==1
	%ni=length(GridDan(idirCCN).Y1)*icediagsALL(idirCCN).i(z0:zend,itdehyd,[226:228])/npess2(idirCCN); %ACC_A * tot no. points = no. icy points Q04-06
    ice=interp1(zz/1000,f*sum(icediagsALL(idirCCN).i(z0:zend,itdehyd,[40:42]),3)/npess2(idirCCN),zzf); %ice value in high CCN case - use min of this and ice sat to moisten towards (i.e. ignore effect of ice above ice sat)
    if ibelsatice==1
		ice=min([ice;sat]); %only use ice if below ice sat
    end
end

if isat==1
    qsat=sat(1); %use first ice sat value as start vapour amount
end
if iqenv==1
    qsat=f*icediagsALL(idir).i(z0,1,37)/npess2(idir); %initial vapour profile value at base - makes this the input base value
end

q(1)=qsat;
q(2:length(zzf))=q(1); %starting values in ppmv - values given to rest of the profile turn out to be irrelevant
%q=repmat(q,[nt 1])';
%q=sat;

if isat==1 & iqenv==0
    [minsat,iminsat]=min(sat);
    q(iminsat:end)=minsat;
end

Nevs=[0:20:240]/10;

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
%     tovflag=0;
%     
%      if tt>=tov
%         tovflag=1; %flag to say it's time to deplete vapour
%         tt=0;
%      end
%      
%     rhoav=(rhof(1:end-1)+rhof(2:end))/2;
    
     
%     for iz=2:length(zzf)
%         
% 
%
    if oldadvect==0
      %   dqdt_as =  ( ww(1:end-1) .* q(1:end-1) - ww(2:end).*q(2:end) ) / dzzf;    %(rho(iz-1) - rho(iz)) * q(iz)/dzz(iz);
          dqdt_as =  ww(2:end) .* ( q(1:end-1) - q(2:end) ) / dzzf;  
     else
         dqdt_as = 0;
     end
    
%         %dqdt_as(iz,t) = dqdt_as(iz,t) - w* (rho(iz-1) - rho(iz)) * q(iz)/dzz(iz);
%         
%  %       if igrad==1 %for gradual removal of vapour - effect spread evenly over time between events (1/N) dq/dt= change/dt = N*change
%            % dqdt_ov(iz,t) =  Nev* prefact(iz)/rhoav/dzz(iz);  %spread based on the sum of dq (dq3) over all times. Done seperately for positive 
%                %dqdt_ov(iz,t) = - Nev*dq_tot(idir).d(iz+z0-1,itdehyd,2);
% %            dqdt_ov(iz,t) = - Nev * ( dq_tot(idir).d(iz+z0-1,itdehyd,2) + nn(idir).n(iz+z0-1,2)/N *(q(iz)-qbase) ); %qbase=q deficit calc'd from e.g.=5ppmv 
%                             %dq/dt=( Nev*dqtot/N + n/N(qenv-5) )    3rd index is 2 as refers to dqtot for <5 ppmv (see Allimp..)
%                             %no need to divide dq_tot by length(Grid.Y1) since was done in Allimp..
%                             
%                 dqdt_ov = - Nev * ( dqtotf(iz) + nntotf(iz)/N *(q(iz)-qbase) ); %qbase=q deficit calc'd from e.g.=5ppmv 
%                 
%                 %       else
% %            dqdt_ov(iz,t)=0;
%             
% %        end
%         
%         
%         q(iz)=q(iz) + ( dqdt_ov ) * dt;
%         q(iz)=max( [0 q(iz)] ); %keep >=0
%         
% %        if tovflag==1 & igrad==0
%                 %q(iz,t)=q(iz,t) + ( AA(iz) - BB(iz)*q(iz,t-1) ); %deplete q if time for overshoot
% %                q(iz,t)=max( [ q(iz,t) - dqz(iz)/sumdqz /rhoav/dzz(iz) 0 ] ); %keep above zero
% %        end
%         
%       
%         
%     end


    
%	dqdt_ov = - Nev * ( dqtotf(2:end)/N + nntotf(2:end)/N/N .*(q(2:end)-qbase) + (n2-N)*q(2:end)/N/N); %qbase=q deficit calc'd from e.g.=5ppmv 
%    dqdt_ov =  Nev* (   n2(2:end)./N .*( qbase*nntotf(2:end)/N - dqtotf(2:end) ) - nntotf(2:end).*n2(2:end)/N/N .* q(2:end)   ); 
%    dqdt_ov =  Nev* ( -q(2:end) + (   pi/4*nntotf(2:end).*(5*nntotf(2:end)-dqtotf(2:end)) + (N*N - pi*nntotf(2:end)/4).*q(2:end) ) );
    dqdt_ov =  Nev*nntotf(2:end)*pi/4/N/N .* ( qbase*nntotf(2:end) - dqtotf(2:end)*N - q(2:end).*nntotf(2:end) );
    
    
    
    if imoisten==1%moisten towards ice sat at a rate of Nmoist events per month
       dqdt_moist =  Nmoist * n2/N* ice(2:end) ; %so not dependant on environmental vapour - dividing by one N less since is mean over all points in icediagsALL
       
       % dqdt_moist =  Nmoist* (   n2/N *( -ice(2:end)/N ) - ni(2:end)*n2/N/N .* q(2:end)   ); %don't need to do this really since are "super-imposing"
                                                        %ice points onto the background field rather than displacing with an amount of tot water
       % dqdt_moist = ( ice(2:end) - q(2:end) ) * Nmoist; %can't really do this as is dependant on q of env - i.e. gets lower as approach ice value
    else
        dqdt_moist=0;
    end
    
    q(2:end)=q(2:end) + ( dqdt_ov + dqdt_moist + dqdt_as ) * dt;
    
    if ibelsat==1
        q=min([q;ones(size(q)).*sat]); %keep q below ice sat
    end
    
    q=max( [q;zeros(size(q))] ); %keep q>=0
    
 if oldadvect==1   
    q(nfi+1:end)=q(1:end-nfi); %nfi is number of points shifted up each timestep
    q(1:nfi)=qsat;         %bring in qbase air from below 
 end
 
 
% qsave(:,t)=q';
    
 %   q(2:100)=4;
    
%break    
end  %t=2:nt

%break

[qmin(n) imin]=min(q(2:end));
hmin(n)=zzf(imin+1);
qq(:,n)=q';

end %n

zzeq=zzf*1000; %for plotting in plotTime..
'done'

infostr=['qbase=' num2str(qbase) ',qsat=' num2str(qsat) ',isatenv=' num2str(isatenv) ',isat=' num2str(isat) ',iqenv=' num2str(iqenv) ...
        ',ibelsat=' num2str(ibelsat) ',ibelsatice=' num2str(ibelsatice) ...
        ',imoisten=' num2str(imoisten) ',Nmoist=' num2str(Nmoist) ',H1=' num2str(H1) ',H2=' num2str(H2) ...
        ',h0=' num2str(h0) ',nclouds=' num2str(nclouds) ',in2EQn1=' num2str(in2EQn1) ',n2=' n2str ',itdehyd=' num2str(itdehyd) ',w=' num2str(w) ...
        ',N1=' num2str(N) ',N2=' num2str(N)];

qstore(istore+1).q=qq;
qstore(istore+1).info=infostr;
istore=istore+1;

