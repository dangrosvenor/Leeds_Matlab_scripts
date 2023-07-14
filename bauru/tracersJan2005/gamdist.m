clear dist

h=16.2;
x=-376.55e3;

minM=10^(-9);
maxM=10^(-6);

isavemulti=0;

comp='lacieLap';
comp='uni';

tag='11thApr06';

iz=findheight(GridDan(idir).Z+620,h*1000);
ix=findheight(GridDan(idir).Y1,x);

%iz=120;
ix=[1:length(GridDan(idir).Y1)];


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

    q=TwoDDan(idir).Q(iz,ix,im);
    n=TwoDDan(idir).Q(iz,ix,in);
    rho=GridDan(idir).RHON(iz);
    [lam,nx0,dist(it).n,dist(it).m,D]=gamlemRow(n,q,rho,lower(lab(1)));
    
%     qstore(it)=q;
%     nstore(it)=n;
%     rhostore(it)=rho;
%     lamstore(it)=lam;
%     nxstore(it)=nx0;
    
    fname=[lab ' mass dist'];
    h1(it).h=figure('name',fname);
    dat=log(10)*D.*dist(it).m.*dist(it).n./rho; %convert from dM/dD to dM/dlogD means *ln(10)*D units: kg(ice)/kg(air) /rho as ndist in #/m/m3
    i=find(dat>minM);
    plot(D(i)*1e6,dat(i));
    
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    %set(gca,'ylim',[minM maxM]);
    xlabel('Diameter (10^{-6} m)'); 
    ylabel([lab ' Mixing Ratio Distribution, dq/dlogD (kg/kg)']); 


    
    savename=[fname ' at x=' num2str(x) ' and h=' num2str(h) 'km'];
     if isavemulti==1
         switch comp
         case 'uni'
             picname=['c:/documents and settings/login/my documents/temp/' savename '_' num2str(idir) '_' tag '.emf'];
         case 'lacieLap'
              picname=['c:/documents and settings/g/my documents/temp/' savename '_' num2str(idir) '_' tag '.emf'];
         end
         print(gcf,picname,'-dmeta');
     end
     
    fname=[lab ' num dist'];
    h2(it).h=figure('name',fname);
    dat=log(10)*D.*dist(it).n./rho; %dN/dlogD units: #/kg
    i=find(dat>1e-3);
    plot(D(i)*1e6,dat(i));
    
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    xlabel('Diameter (10^{-6} m)'); 
    ylabel('Number Distribution dN/dlogD (#/kg)'); 
    
    figs=[figs;h1(it).h;h2(it).h];
    
     savename=[fname ' at x=' num2str(x) ' and h=' num2str(h) 'km'];
     if isavemulti==1
         switch comp
         case 'uni'
             picname=['c:/documents and settings/login/my documents/temp/' savename '_' num2str(idir) '_' tag '.emf'];
         case 'lacieLap'
              picname=['c:/documents and settings/g/my documents/temp/' savename '_' num2str(idir) '_' tag '.emf'];
         end
         print(gcf,picname,'-dmeta');
     end
     
end

tilefigsDAN([3 2],0,[figs]);


[lam,nx0,ns,ms,D]=gamlemRow ( TwoD.Q(100,50,9), TwoD.Q(100,50,4) , GridDan(1).RHON(100),'s');
radar=10*log10( 0.224*sum(D.*ns.*(D*1e3).^6) ); %convert D to mm as Z units in mm6/m3 ns=dN/dD /m3 /m so *D



