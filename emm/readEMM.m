outdir='output';
%outdir='outHM';
outdir='hm2out';
outdir='hm2outMinus50';
outdir='caero1400';


comp='pc';
switch comp
case 'laptop'
	emmdir=['c:/documents and settings/g/my documents/emm_0012/' outdir '/dyn/'];
	emmdir2=['c:/documents and settings/g/my documents/emm_0012/' outdir '/ztime/'];
case 'pc'
    emmdir=['c:/documents and settings/login/my documents/emm_0012/' outdir '/dyn/'];
	emmdir2=['c:/documents and settings/login/my documents/emm_0012/' outdir '/ztime/'];
end

cloudtop=dlmread([emmdir 'cloudtop'],' '); %cldtop,cldtop_max,time
thermalheight=dlmread([emmdir 'thermalheight'],' '); %time, top1, base1, top2,... etc.
thermalwidth=dlmread([emmdir 'thermalwidth'],' ');
W=dlmread([emmdir 'W'],' ');
Wmax=dlmread([emmdir 'Wmax'],' '); %w,alt,time
Wmean_updraft=dlmread([emmdir 'Wmean_updraft'],' ');

a=find(W(:,3)==0.5);
s1=length(a);
b=size(W,1);
s2=b/s1;

wemm=W(:,1);
wemm=reshape(wemm,[s1 s2]); %updraught

temm=W(:,3);
temm=reshape(temm,[s1 s2]); %time grid

hemm=W(:,4);
hemm=reshape(hemm,[s1 s2]); %height grid

Temm=W(:,2);
Temm=reshape(Temm,[s1 s2]); %temp grid

temm=20.25+temm/60;

[temm,Temm,iwczt]=readEMMf(emmdir2,'iwc_up'); %units g/m3   reads in 14 columns with last one as zeroes for some reason (?) - even though only
                                              %13 actually output
[temm,Temm,inczt]=readEMMf(emmdir2,'N_ice_up'); %units per m3
[temm,Temm,naero]=readEMMf(emmdir2,'N_aero_up');% /kg
[temm,Temm,ssat]=readEMMf(emmdir2,'spr_up'); % units percent
[temm,Temm,lwc]=readEMMf(emmdir2,'lwc_up'); %liquid CW (not rain) g/m3
[temm,Temm,ncw]=readEMMf(emmdir2,'N_CW_up'); %#/m3
[temm,Temm,rwc]=readEMMf(emmdir2,'RWC_up'); %rainwater in g/m3
[temm,Temm,nr]=readEMMf(emmdir2,'N_R'); %rainwater concentration /m3



% figure
% 
% pcolor(temm,hemm,wemm);shading interp;colorbar;
% set(gca,'xlim',[GridDan(1).t(1)+3 GridDan(1).t(18)+3]); 
% set(gca,'ylim',[0 22.6969]);

'done'