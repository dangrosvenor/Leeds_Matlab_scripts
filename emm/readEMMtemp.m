clear outdir

starttime=20.15;

    
runs=[9];

outdir{1}='output';
outdir{2}='outHM';
outdir{3}='caero4000_hm2out';
outdir{4}='caero6000_hm2outMinus50';
outdir{5}='caero1400_minus50';
outdir{6}='caero570_2';
outdir{7}='04Jul06';
outdir{8}='22Jul06_nohomog';
outdir{9}='01Aug06_caero4000';

for ii=1:length(runs)
    
    i=runs(ii);
    run_name_emm{ii}=outdir{i};
comp='laptop';
switch comp
case 'laptop'
	emmdir=['c:/documents and settings/g/my documents/emm_0012/' outdir{i} '/dyn/'];
	emmdir2=['c:/documents and settings/g/my documents/emm_0012/' outdir{i} '/ztime/'];
case 'pc'
    emmdir=['c:/documents and settings/login/my documents/emm_0012/' outdir{i} '/dyn/'];
	emmdir2=['c:/documents and settings/login/my documents/emm_0012/' outdir{i} '/ztime/'];
end


cloudtop=dlmread([emmdir 'cloudtop'],' '); %cldtop,cldtop_max,time
thermalheight=dlmread([emmdir 'thermalheight'],' '); %time, top1, base1, top2,... etc.
thermalwidth=dlmread([emmdir 'thermalwidth'],' ');

% ************** just commented this as forgot to copy W *********************
% W=dlmread([emmdir 'W'],' ');
% Wmax=dlmread([emmdir 'Wmax'],' '); %w,alt,time
% Wmean_updraft=dlmread([emmdir 'Wmean_updraft'],' ');
% 
% a=find(W(:,3)==0.5);
% s1=length(a);
% b=size(W,1);
% s2=b/s1;
% 
% wemm=W(:,1);
% wemm=reshape(wemm,[s1 s2]); %updraught
% grid(ii).w=wemm;
% 
% 
% temm2=W(:,3);
% temm2=reshape(temm2,[s1 s2]); %time grid
% temm2=starttime+temm2/60;
% grid(ii).time=temm2;
% 
% Zemm=W(:,4);
% Zemm=reshape(Zemm,[s1 s2]); %height grid
% grid(ii).z=Zemm;
% 
% Temm2=W(:,2);
% Temm2=reshape(Temm2,[s1 s2]); %temp grid
% grid(ii).temp=Temm2;




[temm,Temm,emmdat(ii).iwczt]=readEMMf(emmdir2,'iwc_up',3); %units g/m3   reads in 14 columns with last one as zeroes for some reason (?) - even though only                                              %13 actually output
[temm,Temm,emmdat(ii).inczt]=readEMMf(emmdir2,'N_ice_up',3); %units per m3
[temm,Temm,emmdat(ii).naero]=readEMMf(emmdir2,'N_aero_up',3);% /kg
[temm,Temm,emmdat(ii).ssat]=readEMMf(emmdir2,'spr_up',3); % units percent
[temm,Temm,emmdat(ii).lwc]=readEMMf(emmdir2,'lwc_up',3); %liquid CW (not rain) g/m3
[temm,Temm,emmdat(ii).ncw]=readEMMf(emmdir2,'N_CW_up',3); %#/m3
[temm,Temm,emmdat(ii).rwc]=readEMMf(emmdir2,'RWC_up',3); %rainwater in g/m3
[temm,Temm,emmdat(ii).nr]=readEMMf(emmdir2,'N_R',3); %rainwater concentration /m3
[temm,Temm,emmdat(ii).qisg]=readEMMf(emmdir2,'Q_isg_up',5); %ice+snow+graupel mass g/kg
[temm,Temm,emmdat(ii).qi]=readEMMf(emmdir2,'Q_ice',3); %ice mass g/kg



temm=starttime+temm/60;
% vec(ii).time=temm;
% vec(ii).temp=Temm;
% vec(ii).z=Zemm(:,1);


end

% figure
% 
% pcolor(temm,hemm,wemm);shading interp;colorbar;
% set(gca,'xlim',[GridDan(1).t(1)+3 GridDan(1).t(18)+3]); 
% set(gca,'ylim',[0 22.6969]);

'done'