clear dist

figure

f=1e6*28.97/18;

idir=1;
idir=2;

h=15.5;
h=15.6;
x=-800e3;
ix=findheight(GridDan(idir).Y1,x);

%ix=127;

minM=10^(-9);
maxM=10^(-6);

isavemulti=0;

comp='lacieLap';
comp='uni';

tag='11thApr06';


for idir=1:2

iz=findheight(GridDan(idir).Z+620,h*1000);


tref=repmat(GridDan(idir).THREF(iz,:),[1 length(GridDan(idir).Y1)]);
                       
T=TwoDDan(idir).TH1(iz,:)+tref; %tot potemp

P=TwoDDan(idir).PP(iz,:); %tot P

T=T./(1e5./P).^0.286; %tot temp

rho=P.*28.97e-3/8.3144./T;
            
tot.m=0;
tot.n=0;
tot.dm=0;
tot.dn=0;

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

    q=TwoDDan(idir).Q(iz,:,im);
    n=TwoDDan(idir).Q(iz,:,in);
    q=mean(q);
    n=mean(n);
    
    rho=GridDan(idir).RHON(iz);
    [lam,nx0,dist(it).n,dist(it).m,D]=gamlemRow(n,q,rho,lower(lab(1)));
    
%     qstore(it)=q;
%     nstore(it)=n;
%     rhostore(it)=rho;
%     lamstore(it)=lam;
%     nxstore(it)=nx0;

    tot.m=tot.m+dist(it).m;
    tot.n=tot.n+dist(it).n;
    
    
    tot.dm=tot.dm+dist(it).m.*dist(it).n./repmat(mean(rho),[size(tot.m) 1]);% dq/dD units: kg(ice)/kg(air)/m. /rho as ndist in #/m/m3
	tot.dn=tot.dn+dist(it).n./repmat(mean(rho),[size(tot.m) 1]); %dN/dlogD units: #/kg/m

end



iend=findheight(D,7e-5);   
d=[D(1):D(iend)/500:D(iend)]*1e6;
sum_dm=1e-6*f*interp1(D*1e6,tot.dm,d); 

colstr={'b','r'};
plot(d,sum_dm,colstr{idir});
hold on
end