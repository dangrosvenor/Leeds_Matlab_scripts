%flowing = 0 : debris region (uses w_deb) (if are vertically above thermals but horizontally within width of updraught)
%flowing = 1 : main updraught (uses w_areamean_channel)
%flowing = 2 : outflow region(?) (uses w_out) (scr - if are horizontally outside thermal) 
%flowing = 3 : downdraught?

%function icol=iemm_mphys_cols(dgname,flowing,iflag)
% - use this to find the meaning of the colunmns

comp='lacieLap'; %is called this for consitency with plottimeheight3 even though may not use the Lacie
comp='pc';

source='hd';
%source='lacie';
source='database4';

imphys=0; %additional flag to allow the turning off of all mphysical diag read ins without having to 
            %change the general settings for each run in mphys_diags array

clear outdir

starttime=20.05;
starttime=0;

    
runs=[113 110 111];
%runs=[113 115 116];
runs=113;

outdir{116}='no_homog_aero';
outdir{115}='no_aero_entrain';
outdir{114}='primary_aero_deplete';
outdir{113}='test_cvs_version';
outdir{112}='droplet_test';
outdir{111}='width_2km';
outdir{110}='width_4km';
outdir{109}='noshear_c2000';
outdir{108}='wcrit9_c5000_noshear';
outdir{107}='wcrit9_c5000_2';
outdir{106}='wcrit9_c2000_2';
outdir{105}='wcrit9';
outdir{104}='noshear_dyn_exp';
outdir{103}='wcrit2';
outdir{102}='maxup16_wcrit';
outdir{101}='maxup16_timedecay66';
outdir{100}='maxup16_decay3.0';
outdir{99}='maxup16_cent1.3';
outdir{98}='maxup16_cent1.0';
outdir{97}='norainscr_c2000_wmax30';
outdir{96}='norain_scr_c2000';
outdir{95}='norain_scr';
outdir{94}='mphys_02Aug07_autodiags_noshear_no_liq_scr';
outdir{93}='mphys_02Aug07_autodiags_noshear_norainscr';
outdir{92}='mphys_02Aug07_autodiags_noshear';
outdir{91}='mphys_02Aug07_autodiags';
outdir{90}='mphys_1Aug07_autodiags_0.4';
outdir{89}='mphys_1Aug07_autodiags';
outdir{88}='mphys_31Jul07';
outdir{87}='mphys_diags3';
outdir{86}='c3.0_c5000';
outdir{85}='fifth_c1.5';
outdir{84}='newdyn_fifth';
outdir{83}='fifth_dyn';
outdir{82}='caero5000_0.04';
outdir{81}='nm_newdyn_c5000_0.04';
outdir{80}='new_dyn_c2000_0.04';
outdir{79}='debug';
outdir{78}='nm_new_dyn';
outdir{77}='zero';
outdir{76}='nm_newdyn_test';
outdir{75}='nm2000_0.04';
outdir{74}='nm2000_0.4notbar';
outdir{73}='nm2000_no0.4_noentrain_aero';
outdir{72}='newmex2000_no0.4';
outdir{71}='newmex2000';
outdir{70}='new_mexico';
outdir{69}='noaero_nosec';
outdir{68}='no_aero2';
outdir{67}='meyer_fac_c5000';
outdir{66}='meyer_fac_ztime3';
outdir{65}='glaciation2_c5000';
outdir{64}='glaciation2';
outdir{63}='rain_on_c1000';
outdir{62}='emm_Paul_c1000';
outdir{61}='emm_0013a_Paul';
outdir{60}='old_caero568_test';
outdir{59}='diag_test';
outdir{58}='old_drops_caero1371';
outdir{57}='old_drops_1e-7'; %think is caero=5484
outdir{1}='output';
outdir{2}='outHM';
outdir{3}='caero4000_hm2out';
outdir{4}='caero6000_hm2outMinus50';
outdir{5}='caero1400_minus50';
outdir{6}='caero570_2';
outdir{7}='04Jul06';
outdir{8}='22Jul06_nohomog';
outdir{9}='01Aug06_caero4000';
outdir{10}='dumptest';
outdir{11}='dumptest2';
outdir{12}='ssatdbg';
outdir{13}='zinit1000';
outdir{14}='cent0.9';
outdir{15}='cent2.0';
outdir{16}='new_CB_21Aug06';
outdir{17}='newCB_cent2.0';
outdir{18}='newCB_cent4.0';
outdir{19}='dyntest';
outdir{20}='cent6.0';
outdir{21}='cent4_fe1';
outdir{22}='c4_nofixed';
outdir{23}='c5_nofix_caero1371';
outdir{24}='c5_newent';
outdir{25}='c5_caero5484';
outdir{26}='c1_newQ_1371';
outdir{27}='c5_newQ_1371';
outdir{28}='c1_1371_newQ2';
outdir{29}='c52_1371_newQ2';
outdir{30}='c52_1371_newQ2_nosec';
outdir{31}='c5_caero500_nosec';
outdir{32}='c2_caero500';
outdir{33}='final_c2_caero2000';
outdir{34}='final_c2_caero4000';
outdir{35}='final_c2_caero12000';
outdir{36}='old_droplets_c2_caero2000';
outdir{37}='old_droplets_dis0.1_c2_caero2000';
outdir{38}='old_droplets_caero2000_wide';
outdir{39}='old_drops_dis0.3_caero2000_wide_kesslerON';
outdir{40}='old_drops_dis0.3_caero6000_wide_kesslerON';
outdir{41}='new_drops_caero6000_wide_kesslerON';
outdir{42}='no_sec_new_drops_caero6000_wide_kesslerON';
outdir{43}='no_sec_new_drops_caero20000_wide_kesslerON';
outdir{44}='fth0_new_drops_caero2000_wide_kesslerON';
outdir{45}='fth0_new_drops_caero8000_wide_kesslerON';
outdir{46}='fth0_TIM_NORAIN_caero8000_wide_kesslerON';
outdir{47}='fth0_TIM_NORAIN_3_caero8000_wide_kesslerON';
outdir{48}='sec1st_fth0_TIM_NORAIN_0_caero8000_wide_kesslerON';
outdir{49}='relax5_fth0_TIM_NORAIN_0_caero8000_wide_kesslerON';
outdir{50}='ii_1_minus1_relax5_fth0_TIM_NORAIN_0_caero8000_wide_kesslerON';
outdir{51}='ii_1_plus1_relax5_fth0_TIM_NORAIN_0_caero8000_wide_kesslerON';
outdir{52}='nomerge_relax15_fth0_TIM_NORAIN_0_caero8000_wide_kesslerON';
outdir{53}='zbase_test';
outdir{54}='jjtopminus_1';
outdir{55}='lwc_thresh1e-5';
outdir{56}='threshRAD_1e-10';

gdir=zeros([1 length(outdir)]);
gdir([56 57 58])=1; %ones that are stored on G: drive on PC

%set the runs that have mphys_diags
mphys_diags=zeros([1 length(outdir)]);
%mphys_diags(59)=1;
mphys_diags([87:length(outdir)])=0;

mphys_diags([108 112 113 114])=1;
    





for ii=[1:length(runs)]
%for ii=1:length(runs)
    
    i=runs(ii);
    run_name_emm{ii}=outdir{i};
    



switch comp
case 'lacieLap'
    switch source
    case 'hd'
		rootdir=['c:/documents and settings/g/my documents/emm_0012/' outdir{i}];
    case 'lacie'
        	rootdir=['e:/emm_0012/' outdir{i}];        
    end
case 'pc'
    switch source
    case 'hd'
        if gdir(i)==0;   
			rootdir=['c:/documents and settings/login/my documents/emm_0012/' outdir{i}];
        else
            rootdir=['g:/emm_runs/' outdir{i}];
        end
    case 'lacie'
    	rootdir=['e:/emm_0012/' outdir{i}];       
    case 'database4'
    	rootdir=['y:/' outdir{i}];       
    end
end


    emmdir=[rootdir '/dyn/'];
	emmdir2=[rootdir '/ztime/'];

if imphys==1 & mphys_diags(i)==1
    
    %function icol=iemm_mphys_cols(dgname,flowing,iflag)
    % - use this to find the meaning of the colunmns
    
[mphys(ii).t,mphys(ii).z,mphys(ii).mass]=readEMMf_mphys(emmdir2,'mphys_diags_mass');
[mphys(ii).t,mphys(ii).z,mphys(ii).nums]=readEMMf_mphys(emmdir2,'mphys_diags_nums');

%     mp=dlmread([emmdir2 'mphys_diags'],' ');
%     ncols=size(mp,2)-1;
%     fid=fopen([emmdir2 'mphys_diags'],'rt');
%     mp=fscanf(fid,'%e',[58 Inf]);
%     
%     
%     mp=mp';
%     mp(end,:)=[];
%     
% 	a=find(mp(:,1)==0.5);
% 	s1=length(a);
% 	b=size(mp,1);
% 	s3=size(mp,2);
% 	s2=floor(b/s1)-1;
% 	
% 	mphys=mp(1:s1*s2,:);
% 	mphys=reshape(mphys,[s1 s2 s3]); %re-shape into 3-d array




end
    
    
W=dlmread([emmdir 'W'],' ');
Wmax=dlmread([emmdir 'Wmax'],' '); %w,alt,time
Wmean_updraft=dlmread([emmdir 'Wmean_updraft'],' ');

a=find(W(:,3)==0.5);
s1=length(a);
b=size(W,1);
s2=floor(b/s1);

wemm=W(1:s1*s2,1);
wemm=reshape(wemm,[s1 s2]); %updraught
grid(ii).w=wemm;


temm2=W(1:s1*s2,3);
temm2=reshape(temm2,[s1 s2]); %time grid
temm2=starttime+temm2/60;
grid(ii).time=temm2;

Zemm=W(1:s1*s2,4);
Zemm=reshape(Zemm,[s1 s2]); %height grid
grid(ii).z=Zemm;

Temm2=W(1:s1*s2,2);
Temm2=reshape(Temm2,[s1 s2]); %temp grid
grid(ii).temp=Temm2;

emmdat(ii).W=wemm;
emmdat(ii).T=temm2;
emmdat(ii).Z=Zemm;


rhenv=dlmread([rootdir '/rh_env'],' '); %rh, pressure, temp
emmdat(ii).adiabat=dlmread([rootdir '/adiabat'],' '); %alt=(:,1),pressure=(:,2),temp=(:,3),lwc=(:,4) of adiabat
    
env=dlmread([rootdir '/env'],' ');

cloudtop=dlmread([emmdir 'cloudtop'],' '); %cldtop,cldtop_max,time
thermalheight=dlmread([emmdir 'thermalheight'],' '); %time, top1, base1, top2,... etc.
thermalwidth=dlmread([emmdir 'thermalwidth'],' ');



fprintf(1,'\n\n*********    Finished reading dynamics - ignore future errors if only dynamics required **************\n\n');

[temm,Temm,emmdat(ii).lwc]=readEMMf(emmdir2,'lwc_up',3); %liquid CW (not rain) g/m3

%13 actually output

%temm=starttime+temm/60;

vec(ii).time=unique(temm);
vec(ii).temp=Temm;
vec(ii).z=Zemm(:,1);

[temm,Temm,emmdat(ii).iwczt]=readEMMf(emmdir2,'iwc_up',3); %units g/m3   reads in 14 columns with last one as zeroes for some reason (?) - even though only                                              
[temm,Temm,emmdat(ii).inczt]=readEMMf(emmdir2,'N_ice_up',3); %units per m3
[temm,Temm,emmdat(ii).naero]=readEMMf(emmdir2,'N_aero_up',3);% /kg
[temm,Temm,emmdat(ii).ssat]=readEMMf(emmdir2,'spr_up',3); % units percent
[temm,Temm,emmdat(ii).ncw]=readEMMf(emmdir2,'N_CW_up',3); %#/m3
[temm,Temm,emmdat(ii).rwc]=readEMMf(emmdir2,'RWC_up',3); %rainwater in g/m3
[temm,Temm,emmdat(ii).nr]=readEMMf(emmdir2,'N_R',3); %rainwater concentration /m3
[temm,Temm,emmdat(ii).qisg]=readEMMf(emmdir2,'Q_isg_up',5); %ice+snow+graupel mass g/kg
[temm,Temm,emmdat(ii).qi]=readEMMf(emmdir2,'Q_ice',3); %ice mass g/kg
[temm,Temm,emmdat(ii).qup]=readEMMf(emmdir2,'Q_tracer_up',3); %passive tracer
[temm,Temm,emmdat(ii).lwcscr]=readEMMf(emmdir2,'lwc_scr',3); %passive tracer
[temm,Temm,emmdat(ii).rwcscr]=readEMMf(emmdir2,'rwc_scr',3); %passive tracer
[temm,Temm,emmdat(ii).dcw]=readEMMf(emmdir2,'D_CW',3); %mean liquid diameter (microns)
%[temm,Temm,emmdat(ii).upx]=readEMMf(emmdir2,'updraft_x',5); %position of updraught
[temm,Temm,emmdat(ii).ql]=readEMMf(emmdir2,'QL',3); %mean liquid diameter (microns)
[temm,Temm,emmdat(ii).qr]=readEMMf(emmdir2,'QR',3); %mean liquid diameter (microns)
[temm,Temm,emmdat(ii).iwcscr]=readEMMf(emmdir2,'iwc_scr',3); %units g/m3
[temm,Temm,emmdat(ii).incscr]=readEMMf(emmdir2,'N_ice_scr',3); %units per m3
%[temm,Temm,emmdat(ii).iwcdown]=readEMMf(emmdir2,'iwc_down',3); %units g/m3
[temm,Temm,emmdat(ii).ncwscr]=readEMMf(emmdir2,'N_CW_scr',3); %#/m3

[t_col,Temm,emmdat(ii).precip_col]=readEMMf_col([rootdir '/columnav/'],'precip',1); %#/m3




end

% figure
% 
% pcolor(temm,hemm,wemm);shading interp;colorbar;
% set(gca,'xlim',[GridDan(1).t(1)+3 GridDan(1).t(18)+3]); 
% set(gca,'ylim',[0 22.6969]);

'done'