%calculating effect of differing size distributions on loss from fall speed only
%interactive loss based on a starting distribution - calculate fall speed mass and number loss based on q and n
%adjust q and n based on fall speed loss
%aim is to see if differing distributions make a difference on relevant kind of time scale

clear ice icenc twoD

H=15.25;		
iz=findheight(GridDan(1).Z+620,H*1000);
dz=GridDan(1).Z(iz)-GridDan(1).Z(iz-1);

%find time when have max ice content (take from first run)
[a,b]=max(icediagsALL(1).i(iz,dumprange,42)/npes);
    
idirs=[1 2];
for i=1:length(idirs)
    idir=idirs(i);
    
    %starting values
    ice(1)=icediagsALL(idir).i(iz,b,42)/npes;
    icenc(1)=icediagsALL(idir).i(iz,b,43)/npes;

nsteps=1e4;
dt=10;
tim=[0:dt:dt*nsteps];

    for t=1:length(tim)-1
        twoD.Q(1,1,6)=ice(t);
		twoD.Q(1,1,7)=icenc(1);
		%twoD.Q(1,1,7)=icenc(t); %  NEED TO UNCOMMENT THIS WHEN HAVE NUMBERS        
        Vi=FallSpeed(GridDan(idir).RHON(iz),twoD,'ice',[1 1 1]);
        %VNi=FallSpeed(GridDan(idir),twoD,'ice',[1 1 1]);
        
		dq= -ice(t).*Vi/dz; %flux(kg/m^2)/rho/dz dq in kg/kg/s
	%	dqN= -icenc.*VNi/dz; %flux(#/m^2)/rho/dz dq in kg/kg/s
        
        ice(t+1)=ice(t) + dq*dt;
    %    icenc(t+1)=icenc(t) + dqN*dt;
        
    end

end

ice=ice*f; %convert to ppmv

'finished interactive fall speed'
    