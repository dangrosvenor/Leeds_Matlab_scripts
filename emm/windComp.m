outdir='output';
%outdir='outHM';
outdir='hm2out';

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

temm=19.95+temm/60;

plottimeheightVap3;

xlims=get(gca,'xlim');
ylims=get(gca,'ylim');
clims=get(gca,'clim');
figure
pcolor(temm,hemm,wemm);shading interp;colorbar