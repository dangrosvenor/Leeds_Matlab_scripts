 graph='w';
% graph='naero';
% graph='cloudwater';
% graph='lwc';
% graph='lwc_scr';
% graph='inc';
% graph='incscr';
%  graph='rwc';
%  graph='rwcscr';
%  graph='qr';
% graph='dcw';
%  graph='ncw';
% graph='ncwscr';
%  graph='nr'; 
 graph='iwc';
% graph='iwc_scr';
%graph='qup';
% graph='ssat';

isum=1;

iabc=0;

idat=7; %column of data to plot

add_ground_height=0;

idirsave=idir;

if exist('idirsemm');
	idir=idirsemm(iplot); %emm plot required
end

idirstamp=1;
run_name_emm_select=run_name_emm(idir);

i2d=3;

h1=3.5;
h2=15.5; %in km

if ~exist('vec'); vec.time=0; vec.z=0; end

izmin=findheight(vec(idir).z,h1);
izmax=findheight(vec(idir).z,h2);

if length(izmin)==0; izmin=1; end
if length(izmax)==0; izmax=length(vec(idir).z); end

t1=0;
t2=20.8;
t2=89.5;

if length(vec(idir).time)==1
	it1=1;
    it2=length(emmdat(idir).T(1,:));
else
	it1=findheight(vec(idir).time,t1);
	it2=findheight(vec(idir).time,t2);
    it2=[]; %set to this to show all times
end



if length(it1)==0; it1=1; end
%if length(it2)==0; it2=length(vec(idir).time); t2=vec(idir).time(end); end
if (length(it2)==0)
%    it2=size(emmdat(idir).lwc,2)-2;
    it2=length(vec(idir).time);
    try
        t2=vec(idir).time(it2); 
    catch
        'potential error'
    end
end

%xlabelstr='Time UTC';
xlabelstr='Time (hrs)';

% try
    p = emmdat(idir).adiabat(:,2);
    t = emmdat(idir).adiabat(:,3);
    M = 28*1.67E-27;
	k = 1.38E-23;
	G = 9.81;

    rho=p.*M./k./t;
    
    rho2=interp1(emmdat(idir).adiabat(:,1),rho,vec(idir).z*1000); %interpolate onto the vec grid
    rho=repmat(rho2,[1 length(vec(idir).time)]);
    
% catch
% end


    ixlim=1;    %%%%%%%%%%%%%%%%%%%%%%%%%%%% temporarily setting all xlims  0.85
    xlims=[0.0 1.0];
 %   xlims=[0.0 1.7];
    

switch graph
case 'ssat'
   izovr=1; 
   
   tit(1).tit='Supersaturation (%)';

    pdat(1).p=emmdat(idir).ssat(izmin:izmax,it1:it2,1);
    timesTH(1).t=vec(idir).time(it1:it2);
    zz(1).z=1000*vec(idir).z(izmin:izmax);    
    
    iminovr=0;
	imaxovr=0;
    
    mincovOvr = 0;
    maxcovOvr = 3.2;
    
    clab=1;
    
case 'qup'
   izovr=1; 
   
   tit(1).tit='Tracer MR kg^{-1}';

    pdat(1).p=emmdat(idir).qup(izmin:izmax,it1:it2,1);
    timesTH(1).t=vec(idir).time(it1:it2)/60;
    zz(1).z=1000*vec(idir).z(izmin:izmax);    
    
    iminovr=0;
	imaxovr=0;
    
    mincovOvr = 0;
    maxcovOvr = 3.2;
    maxcovOvr = 2.8;
    
    clab=1;

    
case 'naero'
   izovr=1; 
   
   tit(1).tit='Number of aerosol (cm^{-3})';

    pdat(1).p=emmdat(idir).naero(izmin:izmax,it1:it2,1)/1e6;
    timesTH(1).t=vec(idir).time(it1:it2);
    zz(1).z=1000*vec(idir).z(izmin:izmax);    
    
    iminovr=0;
	imaxovr=0;
    
    mincovOvr = 0;
    maxcovOvr = 3.2;
    maxcovOvr = 2.8;
    
    clab=1;

case 'iwc'
   izovr=1; 
   
   tit(1).tit='Ice Water Content (g kg^{-1})';

    pdat(1).p=emmdat(idir).iwczt(izmin:izmax,it1:it2,1)./rho(izmin:izmax,it1:it2,1);
%    pdat(1).p=emmdat(idir).iwczt(izmin:izmax,it1:it2,1);
    timesTH(1).t=vec(idir).time(it1:it2)/60;
    zz(1).z=1000*vec(idir).z(izmin:izmax);    
    
    iminovr=0;
	imaxovr=-0;
    
%    mincovOvr = 0;
%    maxcovOvr = 14;
%    maxcovOvr = 2.8;
    
    clab=1;
    
case 'iwc_scr'
   izovr=1; 
   
   tit(1).tit='Ice Water Content in SCR(g kg^{-1})';

    pdat(1).p=emmdat(idir).iwcscr(izmin:izmax,it1:it2,1)./rho(izmin:izmax,it1:it2,1);
    timesTH(1).t=vec(idir).time(it1:it2);
    zz(1).z=1000*vec(idir).z(izmin:izmax);    
    
    iminovr=0;
	imaxovr=-0;
    
%    mincovOvr = 0;
%    maxcovOvr = 14;
%    maxcovOvr = 2.8;
    
    clab=1;    
    
case 'ncw'
   izovr=1; 
   
   tit(1).tit='Number of liquid droplets (cm^{-3})';
   tit(1).tit='Number of liquid droplets (kg^{-1}) x 10^6';

    pdat(1).p=emmdat(idir).ncw(izmin:izmax,it1:it2,1)/1e6./rho(izmin:izmax,it1:it2,1);
    
    timesTH(1).t=vec(idir).time(it1:it2)/60;
    zz(1).z=1000*vec(idir).z(izmin:izmax);    
    
    iminovr=0;
	imaxovr=1;
    
    dlogflag=0;
    dlogmin=1e-1;
    
    mincovOvr = 0;
    maxcovOvr = 1200;
%    maxcovOvr = 2.8;
    
    clab=0;
    
case 'ncwscr'
   izovr=1; 
   
   tit(1).tit='Number of liquid droplets in SCR (kg^{-1}) x 10^6';

    pdat(1).p=emmdat(idir).ncwscr(izmin:izmax,it1:it2,1)/1e6./rho(izmin:izmax,it1:it2,1);
    
    timesTH(1).t=vec(idir).time(it1:it2);
    zz(1).z=1000*vec(idir).z(izmin:izmax);    
    
    iminovr=0;
	imaxovr=0;
    
    dlogflag=0;
    dlogmin=1e-1;
    
    mincovOvr = 0;
    maxcovOvr = 3.2;
    maxcovOvr = 2.8;
    
    clab=0;    
    
case 'nr'
   izovr=1; 
   
   tit(1).tit='Number of rain droplets (kg^{-1}) x 10^6';

    pdat(1).p=emmdat(idir).nr(izmin:izmax,it1:it2,1)/1e6./rho(izmin:izmax,it1:it2,1);
    
    timesTH(1).t=vec(idir).time(it1:it2);
    zz(1).z=1000*vec(idir).z(izmin:izmax);    
    
    iminovr=0;
	imaxovr=0;
    
    dlogflag=1;
    dlogmin=1e-1;
    
    mincovOvr = 0;
    maxcovOvr = 3.2;
    maxcovOvr = 2.8;
    
    clab=0;    
    
case 'dcw'
   izovr=1; 
   
   tit(1).tit='Diameter of liquid droplets (microns)';

    pdat(1).p=emmdat(idir).dcw(izmin:izmax,it1:it2,1);
    timesTH(1).t=vec(idir).time(it1:it2);
    zz(1).z=1000*vec(idir).z(izmin:izmax);    
    
    iminovr=0;
	imaxovr=1;
    
    mincovOvr = 0;
    maxcovOvr = 3.2;
    maxcovOvr = 30;
    
    clab=0;
    
case 'qr'
   izovr=1; 
   
   tit(1).tit='RWC (g kg^{-1})';

    pdat(1).p=emmdat(idir).qr(izmin:izmax,it1:it2,1);
    timesTH(1).t=vec(idir).time(it1:it2)/60;
    zz(1).z=1000*vec(idir).z(izmin:izmax);    
    
    iminovr=0;
	imaxovr=0;
    
    mincovOvr = 0;
    maxcovOvr = 3.2;
    maxcovOvr = 4.2;
    
    clab=0;
    
case 'rwc'
   izovr=1; 
   
   tit(1).tit='RWC (g m^{-3})';
   tit(1).tit='RWC (g kg^{-1})';
   
    p = emmdat(idir).adiabat(:,2);
    t = emmdat(idir).adiabat(:,3);
    M = 28*1.67E-27;
	k = 1.38E-23;
	G = 9.81;

    rho=p.*M./k./t;
    
    rho2=interp1(emmdat(idir).adiabat(:,1),rho,vec(idir).z*1000); %interpolate onto the vec grid
    rho=repmat(rho2,[1 length(vec(idir).time)]);
    
    pdat(1).p=emmdat(idir).rwc(izmin:izmax,it1:it2,1)./rho(izmin:izmax,it1:it2);
    timesTH(1).t=vec(idir).time(it1:it2)/60;
    zz(1).z=1000*vec(idir).z(izmin:izmax);    
    
    iminovr=0;
	imaxovr=0;
    
  %   mincovOvr = 1;
%     maxcovOvr = 1;
     maxcovOvr = 15;
    
    clab=0;
        
    

case 'rwcscr'
   izovr=1; 
   
%   tit(1).tit='RWC (g m^{-3})';
  tit(1).tit='RWC (g kg^{-1})';

    p = emmdat(idir).adiabat(:,2);
    t = emmdat(idir).adiabat(:,3);
    M = 28*1.67E-27;
	k = 1.38E-23;
	G = 9.81;

    rho=p.*M./k./t;
    
    rho2=interp1(emmdat(idir).adiabat(:,1),rho,vec(idir).z*1000); %interpolate onto the vec grid
    rho=repmat(rho2,[1 length(vec(idir).time)]);
    
    pdat(1).p=emmdat(idir).rwcscr(izmin:izmax,it1:it2,1)./rho(izmin:izmax,it1:it2);
    timesTH(1).t=vec(idir).time(it1:it2)/60;
    zz(1).z=1000*vec(idir).z(izmin:izmax);    
    
    iminovr=0;
	imaxovr=0;
    
  %   mincovOvr = 1;
%     maxcovOvr = 1;
     maxcovOvr = 15;
    
    clab=0;

    
case 'incscr'
   izovr=1; 
   
%   tit(1).tit='Ice Number Concentration (cm^{-3})';
   tit(1).tit='Ice Number Concentration (kg^{-1} x 10^6)';

    pdat(1).p=emmdat(idir).incscr(izmin:izmax,it1:it2,1)./rho(izmin:izmax,it1:it2,1)/1e6;
%    pdat(1).p=emmdat(idir).incscr(izmin:izmax,it1:it2,1)./rho(izmin:izmax,it1:it2,1)/1e6;    
    timesTH(1).t=vec(idir).time(it1:it2);
    zz(1).z=1000*vec(idir).z(izmin:izmax); 
    
	dlogflag=0;
	dlogmin=1e-3;
        
	iminovr=0;
	imaxovr=1;
    
    mincovOvr = dlog(0,dlogmin);
    maxcovOvr = dlog(1e8,dlogmin);
    
    maxcovOvr = 500;
    clab=0;
    
case 'inc'
   izovr=1; 
   
%   tit(1).tit='Ice Number Concentration (cm^{-3})';
   tit(1).tit='Ice Number Concentration (kg^{-1} x 10^6)';

    pdat(1).p=emmdat(idir).inczt(izmin:izmax,it1:it2,1)./rho(izmin:izmax,it1:it2,1)/1e6;
%    pdat(1).p=emmdat(idir).incscr(izmin:izmax,it1:it2,1)./rho(izmin:izmax,it1:it2,1)/1e6;    
    timesTH(1).t=vec(idir).time(it1:it2)/60;
    zz(1).z=1000*vec(idir).z(izmin:izmax); 
    
	dlogflag=0;
	dlogmin=1e-3;
        
	iminovr=0;
	imaxovr=0;
    
    mincovOvr = dlog(0,dlogmin);
    maxcovOvr = dlog(1e8,dlogmin);
    
    maxcovOvr = 500;
    clab=0;
    
case 'lwc'
   izovr=1; 
   
%   tit(1).tit='LWC (g m^{-3})';
   tit(1).tit='LWC (g kg^{-1})';

    pdat(1).p=emmdat(idir).lwc(izmin:izmax,it1:it2,1)./rho(izmin:izmax,it1:it2);
    timesTH(1).t=vec(idir).time(it1:it2)/60;
    zz(1).z=1000*vec(idir).z(izmin:izmax);    
    
%     iminovr=0;
% 	imaxovr=0;
    
%     mincovOvr = 0;
%     maxcovOvr = 2.6;
 %   maxcovOvr = 2.8;
    
    clab=1;

case 'lwc_scr'
   izovr=1; 
   
   tit(1).tit='LWC (g m^{-3})';

    pdat(1).p=emmdat(idir).lwcscr(izmin:izmax,it1:it2,1);
    timesTH(1).t=vec(idir).time(it1:it2)/60;
    zz(1).z=1000*vec(idir).z(izmin:izmax);    
    
    iminovr=0;
	imaxovr=0;
    
    mincovOvr = 0;
 %   maxcovOvr = 3.2;
 %   maxcovOvr = 2.8;
    
    clab=1;
    
case 'cloudwater'
   izovr=1; 
   
    pdat(1).p=emmdat(idir).iwczt(izmin:izmax,it1:it2,cols(icat));
    timesTH(1).t=vec(idir).time(it1:it2);
    zz(1).z=1000*vec(idir).z;    
    
    clab=0;
    
case 'w'
   izovr=1; 
    tit(1).tit='Updraught (m s^{-1})';
    
%    it1=max( [findheight(emmdat(idir).T(1,:),t1) 1] );
%     if length(it2)==0


%     else    
% 		it2=findheight(emmdat(idir).T(1,:),t2);
%     end

    izmin=findheight(emmdat(idir).Z(:,1),3.5);
    izmax=findheight(emmdat(idir).Z(:,1),16);

    
    timesTH(1).t=emmdat(idir).T(1,:);
    pdat(1).p=emmdat(idir).W(:,:);
    zz(1).z=emmdat(idir).Z(:,1)*1000;
    
    clab=0;
    
    isum=0;
    
%     ixlim=1;
%     xlims=[0.0 1.0];

iylim=1;
iylims=[3.5 15.5];
    
    imaxovr=0;
    maxcovOvr=5;
    
case 'naero'

fact=1e6*28.97/18;
	logflag=1;
	tit(1).tit='N-aero';

    nplots2d=1;
    
    clines=0; %makes black contour lines appear
    clab=0;
    
    i2d=3; %tells it are labelling own x axis
    xlabelstr='Time (sec)';
    
    %minZ=15.8e3;
    %maxZ=17e3;
    ncont=25;
    
    %imaxovr=1;
    maxcovOvr=2.8;
    
    %iminovr=1;
    mincovOvr=2.5;
    
    
   % notsame=1;
    
    izovr=1; %flag to set own z axis
    
    flowing=0;
    
    flow=flowing+1;
    z=squeeze(NN(flow).n(8,:,1));   
		%[z0 zend]=findheight(z,15.8e3,17e3);
        %zz(i).z=z(z0:zend)-620; %620 added later
        zz(i).z=z;
      %  minZ=z(1);
      %  maxZ=z(end);
        
        izmin=1;
        izmax=length(zz(1).z);
        
        pdat(1).p=squeeze(NN(flow).n(idat,:,:));
        imaxovr=[0];
        maxcovOvr=5;
        iminovr=[0];
        mincovOvr=1e8;
        
        timesTH(1).t=squeeze(NN(flow).n(1,1,:))';
        
end

figlab=[tit(1).tit run_name_emm{idir}];
savename=figlab;
idir=idirsave;

if isum==1
    dz=diff(zz(1).z);
    dzz=repmat(dz,[1 size(pdat(1).p,2)]);
    sumpdat=sum( pdat(1).p(2:end,:) .* dzz .* rho(izmin+1:izmax,it1:it2) );
end
