%H=16.25;
%H=17.25;
%H=18;

x=-376.55e3;

clear dist twoD

%idir=3;

npes=npess2(idir);

ice=icediagsALL(idir).i(:,:,42)/npes; %mean ice MR time/height
snow=icediagsALL(idir).i(:,:,40)/npes; %mean snow MR time/height
graupel=icediagsALL(idir).i(:,:,41)/npes; %mean graupel MR time/height

icenc=icediagsALL(idir).i(:,:,43)/npes; %mean ice NC time/height
snownc=icediagsALL(idir).i(:,:,45)/npes; %mean snow NC time/height
graupelnc=icediagsALL(idir).i(:,:,44)/npes; %mean graupel NC time/height

twoD.Q(:,:,6)=ice;
twoD.Q(:,:,7)=icenc;
twoD.Q(:,:,4)=snow;
twoD.Q(:,:,9)=snownc;
twoD.Q(:,:,5)=graupel;
twoD.Q(:,:,8)=graupelnc;


minM=10^(-9);
maxM=10^(-6);

%isavemulti=0;

%comp='lacieLap';
%comp='uni';

%tag='11thApr06';

%H=16.5;
%iz=findheight(GridDan(idir).Z+620,H*1000);

% ix=findheight(GridDan(idir).Y1,x);
% iz=120;


iz=ih;


ix=[1:size(ice,2)];


figs=[];
for it=1:3
    
	switch it
	case 1 %ice
        lab='Ice';
        im=6;
        in=7;
	case 2 %snow
        lab='Snow';
        im=4;
        in=9;
	case 3 %graupel
        lab='Graupel';
        im=5;
        in=8;
	end



    q=twoD.Q(iz,ix,im);
    n=twoD.Q(iz,ix,in);
    rho=GridDan(idir).RHON(iz);
    [lam,nx0,distIce(it).n,distIce(it).m,D]=gamlem(n,q,rho,lower(lab(1)));
    
%    distIce(it).dm=log(10)*repmat(D,[1 size(distIce(it).m,2)]).*distIce(it).m.*distIce(it).n./repmat(rho,[size(distIce(it).m)]);% from dM/dD to dM/dlogD means *ln(10)*D units: kg(ice)/kg(air)/m. /rho as ndist in #/m/m3
	distIce(it).dm=distIce(it).m.*distIce(it).n./repmat(rho,[size(distIce(it).m)]);% dq/dD units: kg(ice)/kg(air)/m. /rho as ndist in #/m/m3

%    distIce(it).dn=log(10)*repmat(D,[1 size(distIce(it).m,2)]).*distIce(it).n./repmat(rho,[size(distIce(it).m)]); %dN/dlogD units: #/kg
    distIce(it).dn=distIce(it).n./repmat(rho,[size(distIce(it).m)]); %dN/dlogD units: #/kg/m

end
    
    
%     fname=[lab ' mass dist'];
%     h1(it).h=figure('name',fname);
 %   dat=log(10)*D.*dist(it).m.*dist(it).n./rho; %convert from dM/dD to dM/dlogD means *ln(10)*D units: kg(ice)/kg(air) /rho as ndist in #/m/m3
%     i=find(dat>minM);
%     plot(D(i)*1e6,dat(i));
    
%     set(gca,'xscale','log');
%     set(gca,'yscale','log');
%     %set(gca,'ylim',[minM maxM]);
%     xlabel('Diameter (10^{-6} m)'); 
%     ylabel([lab ' Mixing Ratio Distribution, dq/dlogD (kg/kg)']); 


    
%     savename=[fname ' at x=' num2str(x) ' and h=' num2str(h) 'km'];
%      if isavemulti==1
%          switch comp
%          case 'uni'
%              picname=['c:/documents and settings/login/my documents/temp/' savename '_' num2str(idir) '_' tag '.emf'];
%          case 'lacieLap'
%               picname=['c:/documents and settings/g/my documents/temp/' savename '_' num2str(idir) '_' tag '.emf'];
%          end
%          print(gcf,picname,'-dmeta');
%      end
%      
%     fname=[lab ' num dist'];
%     h2(it).h=figure('name',fname);
%     dat=log(10)*D.*dist(it).n./rho; %dN/dlogD units: #/kg
%     i=find(dat>1e-3);
%     plot(D(i)*1e6,dat(i));
%     
%     set(gca,'xscale','log');
%     set(gca,'yscale','log');
%     xlabel('Diameter (10^{-6} m)'); 
%     ylabel('Number Distribution dN/dlogD (#/kg)'); 
%     
%     figs=[figs;h1(it).h;h2(it).h];
%     
%      savename=[fname ' at x=' num2str(x) ' and h=' num2str(h) 'km'];
%      if isavemulti==1
%          switch comp
%          case 'uni'
%              picname=['c:/documents and settings/login/my documents/temp/' savename '_' num2str(idir) '_' tag '.emf'];
%          case 'lacieLap'
%               picname=['c:/documents and settings/g/my documents/temp/' savename '_' num2str(idir) '_' tag '.emf'];
%          end
%          print(gcf,picname,'-dmeta');
%      end
%      
% end
% 
% tilefigsDAN([3 2],0,[figs]);



