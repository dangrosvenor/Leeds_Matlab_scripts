%quick calculation to work out mean reduction in MR over 1.2km region (i.e. what stratospheric reduction would
%be if there was one event per month with 1.2 km/month ascent rate

f=1e6*28.97/18;

[iz1,iz2]=findheight(GridDan(idir).Z/1000+0.62,15.8,17);
[iz1,iz2]=findheight(GridDan(idir).Z/1000+0.62,15.8,30);
iz2=250;

tot=f*sum(icediagsALL(1).i(iz1:iz2,:,[37:42]),3)/npess2(idir);
z=GridDan(1).Z(iz1:iz2)+620;
dz=z(2:end)-z(1:end-1); 
rho=GridDan(idir).RHO(iz1:iz2);

rho=repmat(rho,[1 size(icediagsALL(idir).i,2)]);
dz=repmat(dz,[1 size(icediagsALL(idir).i,2)]);

vapmass=sum( rho(2:end,:) .*dz .* tot(2:end,:) );
airmass=sum( rho(2:end,:) .*dz );

meanq=vapmass./airmass;
    