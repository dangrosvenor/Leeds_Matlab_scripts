it=2;

ih=9;
ih2=17;

sq=size(ecmwf(1).q);

meecm_simple=mean( squeeze(ecmwf(1).q(it,ih:ih2,:,:)) );

p=repmat(flipud(ecmwf(1).p),[1 25 27]);
t=squeeze(ecmwf(1).t(it,:,:,:));
rho=p*100.*28.97e-3/8.3144./t;
dz=repmat(diff(z)',[1 sq(3) sq(4)]);

meecm=sum( dz(ih:ih2,:,:).*rho(ih:ih2,:,:).*squeeze(ecmwf(1).q(it,ih:ih2,:,:)) );
air=sum(dz(ih:ih2,:,:).*rho(ih:ih2,:,:));

meecm2=meecm./air;

meecmT_simple=mean( squeeze(ecmwf(1).t(it,ih:ih2,:,:)) );
meecmT=sum( dz(ih:ih2,:,:).*rho(ih:ih2,:,:).*squeeze(ecmwf(1).t(it,ih:ih2,:,:)) ) ./ air;

ihh=1;
ihh2=29;

melem_simple=mean(icediagsALL(1).i(ihh:ihh2,:,37)/npess2(1));

dz2=repmat(diff(GridDan(1).Z),[1 size(icediagsALL(idir).i,2)]);
rho2=repmat(GridDan(1).RHON,[1 size(icediagsALL(idir).i,2)]);
melem=sum( rho2(ihh:ihh2,:).*dz2(ihh:ihh2,:).*icediagsALL(idir).i(ihh:ihh2,:,37)/npess2(1) );
air2=sum(dz2(ihh:ihh2,:).*rho2(ihh:ihh2,:));

melem2=melem./air2;


p=repmat(GridDan(1).PREFN,[1 size(icediagsALL(idir).i,2)]);
th=icediagsALL(1).i(:,:,246)/npess2(1);
T=th./(1000e2./p).^0.286;
melemT_simple=mean( T(ihh:ihh2,:) );

melemT=sum( rho2(ihh:ihh2,:).*dz2(ihh:ihh2,:).*T(ihh:ihh2,:) ) ./ air2;




latB=findheight(ecmwf(1).lat,-22.36);
lonB=findheight(ecmwf(1).lon-360,-49.03);

qb=meecm_simple(1,latB,lonB);

%diff for qv of 2.25 g/kg occurs at lat=15, lon=11 = -20, -52
%diff for qv of 2.25 g/kg occurs at lat=15, lon=11 = -20, -52
%diff for qv of 2.25 g/kg occurs at lat=15, lon=8 = -20, -55



