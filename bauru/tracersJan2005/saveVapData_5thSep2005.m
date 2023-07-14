%%%%%%%%%%% Don't put numbers at the end of the filename or matlab will need a .mat at the end ***********************

comp='uni';
%comp='lacieLap';

ndir=99;

switch comp
	case 'uni'
        topdir(1).dir='c:/cygwin/home/user/runs/';
        topdir(2).dir='c:/cygwin/home/user/runs/';
        topdir(3).dir='c:/cygwin/home/user/runs/';
        topdir(4).dir='g:/runs/';
        topdir(5).dir='h:/';
        topdir(6).dir='g:/runs/';
        topdir(7).dir='h:/';
        topdir(8).dir='g:/runs/';
        topdir(9).dir='';
        topdir(10).dir='';
        topdir(11).dir='';
        topdir(12).dir='';
        topdir(13).dir='h:/';
        topdir(14).dir='h:/';
        topdir(15).dir='g:/runs/';
        topdir(16).dir='c:/cygwin/home/user/runs/';
        topdir(17).dir='c:/cygwin/home/user/runs/';
        topdir(18).dir='c:/cygwin/home/user/runs/';
        topdir(19).dir='c:/cygwin/home/user/runs/';
        topdir(20).dir='c:/cygwin/home/user/runs/';
        topdir(21).dir='c:/cygwin/home/user/runs/';
        topdir(22).dir='c:/cygwin/home/user/runs/';
        topdir(23).dir='c:/cygwin/home/user/runs/';
        topdir(24).dir='c:/cygwin/home/user/runs/';
        topdir(25).dir='h:/';
        topdir(26).dir='h:/';
        topdir(27).dir='h:/';
        topdir(28).dir='c:/cygwin/home/user/runs/';
        topdir(29).dir='c:/cygwin/home/user/runs/';
        topdir(30).dir='g:/runs/';
        topdir(31).dir='g:/runs/';
        topdir(32).dir='g:/runs/';
        topdir(33).dir='g:/runs/';
        topdir(34).dir='g:/runs/';  
        topdir(35).dir='g:/runs/';
        topdir(36).dir='g:/runs/';
        topdir(37).dir='g:/runs/';  
        topdir(39).dir='g:/runs/';  
        topdir(40).dir='g:/runs/';  
        topdir(41).dir='g:/runs/';  
        topdir(42).dir='g:/runs/';  
        topdir(43).dir='g:/runs/';  
        topdir(44).dir='g:/runs/';
        topdir(45).dir='g:/arm_radar/echo_tops/';
        topdir(46).dir='g:/runs/';
        topdir(47).dir='g:/runs/';
        topdir(48).dir='g:/runs/';
        topdir(49).dir='g:/runs/';
        topdir(50).dir='g:/runs/';
        topdir(51).dir='g:/runs/';
        topdir(52).dir='g:/runs/';
        topdir(53).dir='g:/runs/';
        topdir(54).dir='g:/runs/';
        topdir(55).dir='g:/runs/';
        topdir(56).dir='g:/runs/';
        topdir(57).dir='g:/runs/';
        topdir(58).dir='g:/runs/';
        topdir(59).dir='g:/runs/';
        topdir(60).dir='g:/runs/';
        topdir(61).dir='g:/runs/';
        topdir(62).dir='g:/runs/';
        topdir(63).dir='g:/runs/';
        topdir(64).dir='g:/runs/';
        topdir(65).dir='g:/runs/';
        topdir(66).dir='g:/runs/';
        topdir(67).dir='g:/runs/';
        topdir(68).dir='g:/runs/';
        topdir(69).dir='g:/runs/';
        topdir(70).dir='g:/runs/';
        topdir(71).dir='g:/runs/';
        topdir(72).dir='g:/runs/';
        
        topdir(73).dir='g:/runs/';        
%        topdir(73).dir='z:/';        
        
        topdir(74).dir='g:/runs/';
        topdir(75).dir='g:/runs/';
        topdir(76).dir='z:/';
        topdir(77).dir='g:/runs/';
        topdir(78).dir='z:/';
        topdir(79).dir='z:/'; 
        topdir(80).dir='z:/';
        topdir(81).dir='z:/'; 

        for itopdir=82:140
            topdir(itopdir).dir='z:/'; 
        end        
        
        for irun=[112 113 118 141:200];
            topdir(irun).dir='y:/LEM/';
        end
        
        
    case 'lacieLap'
        for i=1:ndir
            topdir(i).dir='e:/les/';
        end
        
end

direc(1).dir=[topdir(1).dir 'dmidamp_2/'];

direc(2).dir='c:/cygwin/home/user/runs/damp_ccn960/';
direc(3).dir='c:/cygwin/home/user/runs/damp_inx10/';

direc(4).dir='g:/runs/500mres/';

direc(5).dir='h:/5ppmv_qnums2/';

%direc(1).dir='c:/cygwin/home/user/runs/5ppmv/';
direc(6).dir='g:/runs/5ppmv_ccn960_2/';
direc(7).dir='g:/runs/5ppmv/';

direc(8).dir='h:/5ppmv_1000km_3/';

direc(9).dir='g:/runs/4ppmv_th26.8_2/';
direc(10).dir='g:/runs/5ppmv_th26.8_2/';
direc(11).dir='g:/runs/5ppmv_indiv10/';

direc(12).dir=[topdir(12).dir '4ppmv/']; %

direc(13).dir='h:/lowE/';

direc(14).dir=[topdir(14).dir '250mres_1000km/']; % **zdmp=22.858km
direc(15).dir=[topdir(15).dir '3dNewton/']; 
direc(16).dir=[topdir(16).dir '500mNewt/']; 
direc(17).dir=[topdir(17).dir '1kmNewt/']; 

direc(18).dir=[topdir(18).dir '250m13.4Newt/']; %1000km

direc(19).dir=[topdir(19).dir '13.02.12utc/']; 
direc(20).dir=[topdir(19).dir '13.02_2/']; 

direc(20).dir=[topdir(20).dir '13.02_2/']; %as above but with warmpool at x=10km so all of it will be within domain - both done on PC
direc(21).dir=[topdir(21).dir '13.02_inx10/'];

direc(22).dir=[topdir(22).dir '13.02_ccn960/'];

direc(23).dir=[topdir(23).dir 'ccn960Newt/']; 
direc(24).dir=[topdir(24).dir 'inx10newt/'];

direc(25).dir=[topdir(25).dir 'MilesCity/'];

direc(26).dir=[topdir(26).dir 'Hector_Paul/'];
direc(27).dir=[topdir(27).dir 'Horace_z150_2/'];
direc(28).dir=[topdir(28).dir '1km_norad/'];
direc(29).dir=[topdir(29).dir '1kmccn960/'];
direc(30).dir=[topdir(30).dir 'ccn960_2000km/'];

direc(32).dir=[topdir(32).dir 'damp_10secs_8km/']; 
direc(33).dir=[topdir(33).dir 'iceacc_1000km/']; 

direc(34).dir=[topdir(34).dir 'MilesCity13.4/'];
direc(35).dir=[topdir(35).dir 'Miles13.4_xr2000/'];
direc(36).dir=[topdir(36).dir 'Miles13.4_dry_bot/'];

direc(37).dir=[topdir(37).dir '2000km/'];

direc(39).dir=[topdir(39).dir '4000km/'];

direc(40).dir=[topdir(40).dir 'milescity13.4_ccn960/'];

direc(41).dir=[topdir(41).dir '2000km_20mins/'];
direc(42).dir=[topdir(42).dir 'hightop/'];
direc(43).dir=[topdir(43).dir '2000km_ccn960_2/'];
direc(44).dir=[topdir(44).dir '2000km_xr3500/'];
direc(45).dir=[topdir(45).dir '20051116/'];
direc(46).dir=[topdir(46).dir '2000km_13.4th/'];
direc(47).dir=[topdir(47).dir 'z_150/'];
direc(48).dir=[topdir(48).dir '2km_2000km/'];
direc(49).dir=[topdir(49).dir '2km/'];
direc(50).dir=[topdir(50).dir '4000km_2/'];
direc(51).dir=[topdir(51).dir '3d_sm30/'];
direc(52).dir=[topdir(52).dir '2000km_6.7th/'];
direc(53).dir=[topdir(53).dir '2000km_6.7th_noqv/'];
direc(54).dir=[topdir(54).dir '1.25km/'];
direc(55).dir=[topdir(55).dir 'Redel_20mins/'];
direc(56).dir=[topdir(56).dir '6.7th_1.67qv/'];
direc(57).dir=[topdir(57).dir 'highRH/'];
direc(58).dir=[topdir(58).dir 'nodry/'];
direc(59).dir=[topdir(59).dir '3d2km_30red/'];
direc(60).dir=[topdir(60).dir 'tgraup2/'];
direc(61).dir=[topdir(61).dir '4000km_3.5km/'];
direc(62).dir=[topdir(62).dir 'convap_allmoist/'];
direc(63).dir=[topdir(63).dir '6.7th_20mins/'];
direc(64).dir=[topdir(64).dir 'nodry_highRH/'];
direc(65).dir=[topdir(65).dir 'moist_noqv/'];
direc(66).dir=[topdir(66).dir '3d2km_30/'];
direc(67).dir=[topdir(67).dir 'highT_test/'];
direc(68).dir=[topdir(68).dir 'highT_1.5/'];
direc(69).dir=[topdir(69).dir 'highT1.5_6.7th/'];
direc(70).dir=[topdir(70).dir 'noshear/'];
direc(71).dir=[topdir(71).dir '2km_20mins/'];
direc(72).dir=[topdir(72).dir '6.7th_20mins/'];

direc(73).dir=[topdir(73).dir '3d2km_30_20mins/'];

direc(74).dir=[topdir(74).dir 'ccn960_20m/'];
direc(75).dir=[topdir(75).dir '3d2km_30_20mins/'];
direc(76).dir=[topdir(76).dir '3d2km_6.7th_20mins/'];
direc(77).dir=[topdir(77).dir '3d2km_30red/']; %6.7th and 1.67qv
direc(78).dir=[topdir(78).dir '3d2km_30red/']; %6.7th and 1.67qv
direc(79).dir=[topdir(length(direc)).dir '3d_75km/']; %6.7th and 1.67qv
direc(80).dir=[topdir(length(direc)).dir 'NM_lowf_1e-3_20km/']; 
direc(81).dir=[topdir(length(direc)).dir 'NM_lowf_1e-3_15km/']; 

irun=105;  direc(irun).dir=[topdir(irun).dir 'NM_runs/2e-3_ccn/NM_2e-3_240ccn/']; 
irun=106;  direc(irun).dir=[topdir(irun).dir 'NM_runs/2e-3_ccn/NM_2e-3_480ccn/']; 
irun=107;  direc(irun).dir=[topdir(irun).dir 'NM_runs/2e-3_ccn/NM_2e-3_720ccn/']; 
irun=108;  direc(irun).dir=[topdir(irun).dir 'NM_runs/2e-3_ccn/NM_2e-3_960ccn/']; 


irun=110;  direc(irun).dir=[topdir(irun).dir '2000km_26.8th/']; 
irun=111;  direc(irun).dir=[topdir(irun).dir '2000km_33.5th/']; 

irun=112;  direc(irun).dir=[topdir(irun).dir '150km_2.23th/']; 
irun=113;  direc(irun).dir=[topdir(irun).dir '150km_4.46th/']; 

irun=114;  direc(irun).dir=[topdir(irun).dir '2000km_1e-3_1hour/']; 
irun=115;  direc(irun).dir=[topdir(irun).dir '2000km_2e-3_1hour/']; 
irun=116;  direc(irun).dir=[topdir(irun).dir '2000km_3e-3_1hour/']; 

irun=118;  direc(irun).dir=[topdir(irun).dir '150km_trip20.1/']; 

irun=127;  direc(irun).dir=[topdir(irun).dir 'texas_superad/']; 
irun=128;  direc(irun).dir=[topdir(irun).dir 'texas_superad_1e-3/']; 
irun=129;  direc(irun).dir=[topdir(irun).dir 'newmexico_superad/']; 

irun=130;  direc(irun).dir=[topdir(irun).dir 'newmexico_superad_1e-3/']; 
irun=131;  direc(irun).dir=[topdir(irun).dir 'toga_1e-3_1hr/']; 

irun=133;  direc(irun).dir=[topdir(irun).dir '150km_ccn960_non_trip/']; 

irun=141;  direc(irun).dir=[topdir(irun).dir '75km_bbigg0/']; %
irun=142;  direc(irun).dir=[topdir(irun).dir '75km_orig_sounding/']; %

irun=148;  direc(irun).dir=[topdir(irun).dir '75km_orig_hotpool/']; %

irun=156;  direc(irun).dir=[topdir(irun).dir 'nm_supp_2deg_hightop/']; %

irun=167;  direc(irun).dir=[topdir(irun).dir '150km_ccn960_homog_only/']; %


endt={'65' '88' '65' '88' '58' '88' '55' '60' '' '' '' '88' '' '' '' '' '' ''};
%endtALL={'65' '65' '65' '88' '58' '88' '55' '60' '' '' '' '88' '' '' '' '' '' ''};
endtWmax={'88' '88' '88' '88' '88' 'xx' '55' '60' '' '' '' '88' '' '' '' '' '' ''};

%select things to save

saveflag=zeros([1 999]);

%saveflag([4 5 6 7 8 9 12 17 29])=1; % usual items to save for older cases with no num_rates
%saveflag([4 5 6 7 8 11 12 17 29 39])=1; % usual items to save for new cases with num_rates%

%saveflag([5 6 7 8 11 39])=1; % usual items for limited cases


saveflag(1)=1; %for grid
%saveflag([4:8 11 12 17 29 34])=1; %ususal suite of diags

%saveflag([52 62 8 11 262])=1; %ususal suite of diags **for 3D** run no tot
%saveflag([5 6 8 11 26])=1; %ususal suite of diags **for 3D** run with tot
%saveflag([8 11])=1; %ususal suite of diags **for 3D** run no tot

%saveflag([8 11])=1; %icediagsALL and icediags_nums
%saveflag([8])=1; %icediagsALL'y'



%saveflag(28)=1; %signTpert_lnb
%saveflag([4])=1; %maxW
%saveflag([6])=1; %vap_prctiles
%saveflag(11)=1;
%saveflag([15])=1; %radar timH
%saveflag(16)=1;
%saveflag(21)=1; %lnbbins for vap and tot < 5ppmv
%saveflag(24)=1; %temperature diags
%saveflag(25)=1; %rho5ppmv etc.
%saveflag(26)=1; %vapdist.
%saveflag(30)=1; %icediagsALL for just dump 62 (in order to get ACC_A for eq. model)
%saveflag(31)=1; %TRACERacc3 etc.
%saveflag(32)=1; %INCacc_alltim etc.
%saveflag(33)=1; %sf4 soundings
%saveflag(35)=1; %radar 10 dBz number of points   n10dbz
%saveflag(36)=1; %ARM Darwin radar echo top stats
%saveflag(36)=1; %ARM Darwin radar echo top stats
%saveflag(37)=1; %temperature perturbations of full 2d field
%saveflag(38)=1; %cape fields for grid profiles
%saveflag(39)=1; %pos and neg q' contribution'y's in up and downdraughts for vapour and total water
%saveflag(40)=1; %temperature perturbations of full 2d field
%saveflag(42)=1; %horizontal slices at one height of temp, pressure and vapour for si calc
%saveflag(43)=1; %horizontal slices at one height of iceMR, snowMR, snowNC, snowMR
%saveflag(44)=1; %no. points at each height above threshold lwc
%saveflag(45)=1; %tot ice or water
%saveflag(46)=1; %temp change for conv and non-conv regions. 
%saveflag(47)=1; %'tpertTimH_full','tpert_max' - temperature perts from median for all of heating area and max value at each height

saveselect=226; %directory to save to


        

for jdir=1:length(exdirSTORE)   %length(direc)
    exdir=exdirSTORE(jdir).dat
  
        
    j=saveselect(jdir);
    if j>12
        endt{j}='';
        endtALL{j}='';
        endtWmax{j}='';
    end
    
    %exdir=[direc(j).dir 'results/diags/'];
    
    yn=input(['Proceed with save in directory ' exdir ' ??? (''y'' for yes) : ']);
    if strcmp(yn,'y')
        fprintf(1,'\nSaving...');

    else
        fprintf(1,'\nFiles NOT SAVED');
        break
    end
    
    
    if saveflag(1)==1
		exname=[exdir 'gridDan.mat'];
		gridDan=GridDan(jdir);
		%gridDan.t=time;
		save(exname,'gridDan');
	end
		
	if saveflag(2)==1
		exname=[exdir 'icediag4_5thSep2005'];
		icediags4=icediag4(jdir);
		save(exname,'icediags4');
	end
	
	if saveflag(3)==1
		exname=[exdir 'icediag_5thSep2005'];
		icediags=icediag(jdir);
		save(exname,'icediags');
	end
	
	if saveflag(4)==1
        %exname=[exdir 'maxW_1-' endt{j}];
        exname=[exdir 'maxW.mat'];
        
        maxW=MaxW(jdir);
        minW=MinW(jdir);

        save(exname,'maxW','minW');
        
	end    
    
    if saveflag(5)==1
%		exname=[exdir 'dq_dehyd_1-' endt{j}];
		exname=[exdir 'dq_dehyd.mat'];
		dq_tots=dq_tot(jdir);
        dq_vap=dq_vaps(jdir);
        save(exname,'dq_vap','dq_tots');
%        save(exname,'dq_vap');
        
	end
    
    if saveflag(52)==1
%		exname=[exdir 'dq_dehyd_1-' endt{j}];
		exname=[exdir 'dq_dehyd.mat'];
%		dq_tots=dq_tot(jdir);
        dq_vap=dq_vaps(jdir);
%        save(exname,'dq_vap','dq_tots');
        save(exname,'dq_vap');
        
	end
    
    if saveflag(6)==1
%		exname=[exdir 'vap+tot_prcs_1-' endt{j}];
		exname=[exdir 'vap+tot_prcs.mat'];
		vap_prc=vap_prctiles(jdir);
        tot_prc=tot_prctiles(jdir);
		save(exname,'tot_prc','vap_prc');
%		save(exname,'vap_prc');
        
	end
    
    if saveflag(62)==1
%		exname=[exdir 'vap+tot_prcs_1-' endt{j}];
		exname=[exdir 'vap+tot_prcs.mat'];
		vap_prc=vap_prctiles(jdir);
       % tot_prc=tot_prctiles(jdir);
%		save(exname,'tot_prc','vap_prc');
		save(exname,'vap_prc');
        
	end
    
    if saveflag(7)==1
%		exname=[exdir 'nn_1-' endt{j}];
		exname=[exdir 'nn.mat'];
        n=nn(jdir);
        n2=nn2(jdir);
		save(exname,'n','n2');
%		save(exname,'n2');
        
	end
    
    if saveflag(8)==1
%		exname=[exdir 'icediagsALL_1-' endt{j}];
		exname=[exdir 'icediagsALL.mat'];
        icediagALL=icediagsALL(jdir);
        save(exname,'icediagALL');
        
	end
      
    if saveflag(9)==1
%		exname=[exdir 'icediag_1-' endt{j}];
		exname=[exdir 'icediag'];
        icediags=icediag(jdir);
		save(exname,'icediags');
	end
    
    if saveflag(10)==1
        clear flux
		exname=[exdir 'fluxes_1-65'];
		flux(1).f=icediagALL(jdir).i(:,:,136:138);
		save(exname,'flux');
	end
    
    if saveflag(11)==1
%		exname=[exdir 'icediag_nums_1-' endt{j}];
		exname=[exdir 'icediag_nums.mat'];
        icediags_nums=icediag(jdir);
		save(exname,'icediags_nums');
	end
    
    if saveflag(12)==1
        %exname=[exdir 'maxW_1-' endt{j}];
        clear wd wb
        exname=[exdir 'W_prc+dist.mat'];
        w=w_prctiles(jdir);
        wd=wdist(jdir);
        wb=wbb(jdir);
        
        save(exname,'w','wb','wd');
        
    end
        
    if saveflag(13)==1
        %exname=[exdir 'maxW_1-' endt{j}];
        clear wd wb
        
        Wind=wind(jdir);
        exname=[exdir 'wind_timH_1-22'];
        save(exname,'Wind');
        
        Vap=vap(jdir);
        exname=[exdir 'vap_timH_1-22'];
        save(exname,'Vap');
        
        Pressure=pressure(jdir);
        exname=[exdir 'pressure_timH_1-22'];
        save(exname,'Pressure');
        
        Potemp=potemp(jdir);
        exname=[exdir 'potemp_timH_1-22'];
        save(exname,'Potemp');
        
        IceMR=icemr(jdir);
        exname=[exdir 'iceMR_timH_1-22'];
        save(exname,'IceMR');
        
        IceNC=icenc(jdir);
        exname=[exdir 'iceNC_timH_1-22'];
        save(exname,'IceNC');
        
        SnowMR=snowmr(jdir);
        exname=[exdir 'snowMR_timH_1-22'];
        save(exname,'SnowMR');
        
        SnowNC=snownc(jdir);
        exname=[exdir 'snowNC_timH_1-22'];
        save(exname,'SnowNC');
        
        GraupelMR=graupelmr(jdir);
        exname=[exdir 'graupelMR_timH_1-22'];
        save(exname,'GraupelMR');
        
        GraupelNC=graupelnc(jdir);
        exname=[exdir 'graupelNC_timH_1-22'];
        save(exname,'GraupelNC');
        
	end
    
    if saveflag(15)==1
        clear zm
        exname=[exdir 'radar_max_timH.mat'];
        zm=zmax(jdir);        
        save(exname,'zm');
    end
    
    if saveflag(16)==1
        
        dqnonSAVE=dqnon(jdir);
        nnnonSAVE=nnnon(jdir);
        
        dqpoSAVE=dqpo(jdir);
        nnpoSAVE=nnpo(jdir);
                    
                    
        exname=[exdir 'dq_potemp_medacc0.05.mat'];
        
        save(exname,'dqnonSAVE','nnnonSAVE','dqpoSAVE','nnpoSAVE');
    end
    
    if saveflag(17)==1
        
        simaxSAVE=simaxTimH(jdir);
        %siminSAVE=siminTimH(jdir);
        %simeanSAVE=simean(jdir);
        
        exname=[exdir 'simaxTimH2.mat'];
        
        %save(exname,'simaxSAVE','simeanSAVE','siminSAVE');
        save(exname,'simaxSAVE');
    end    
    
    if saveflag(18)==1
        
        tpertSAVE=tpertTimH(jdir);
        
        exname=[exdir 'Tpert.mat'];
        
        save(exname,'tpertSAVE');
        

    end
    
    if saveflag(19)==1
        
        rhopertmaxSAVE2=rhopertTimHmax2(jdir);
        rhopertminSAVE2=rhopertTimHmin2(jdir);
        rhopertSAVE2=rhopertTimH2(jdir);
        
        rhopertmaxSAVE=rhopertTimHmax(jdir);
        rhopertminSAVE=rhopertTimHmin(jdir);
        rhopertSAVE=rhopertTimH(jdir);

       % tpertSAVE=tpertTimH(jdir);
        
        exname=[exdir 'rhopertALL2.mat'];
        
%        save(exname,'rhopertSAVE','tpertSAVE');
        
         save(exname,'rhopertSAVE','rhopertmaxSAVE','rhopertminSAVE','rhopertSAVE2','rhopertmaxSAVE2','rhopertminSAVE2')
     end
         
    if saveflag(20)==1
        
        me_lnb_abvSAVE=meanlnb_abv(jdir);
        me_lnb_belSAVE=meanlnb_bel(jdir);
        
        minlnbSAVE=minlnb(jdir);
        maxlnbSAVE=maxlnb(jdir);
        
        lnbbins_posSAVE=lnbbins_pos(jdir);
        lnbbins_negSAVE=lnbbins_neg(jdir);
        binsSAVE=bins(jdir).b;
        
        
        
        exname=[exdir 'LNB.mat'];
        
        save(exname,'me_lnb_abvSAVE','me_lnb_belSAVE','minlnbSAVE','maxlnbSAVE','lnbbins_posSAVE','lnbbins_negSAVE')
         

    end
    
    if saveflag(21)==1
        
        me_lnb_abv_vapSAVE=meanlnb_abv_vap(jdir);
        me_lnb_bel_vapSAVE=meanlnb_bel_vap(jdir);
        
        minlnb_vapSAVE=minlnb_vap(jdir);
        maxlnb_vapSAVE=maxlnb_vap(jdir);
        
        lnbbins_pos_vapSAVE=lnbbins_pos_vap(jdir);
        lnbbins_neg_vapSAVE=lnbbins_neg_vap(jdir);
        
        me_lnb_abv_totSAVE=meanlnb_abv_tot(jdir);
        me_lnb_bel_totSAVE=meanlnb_bel_tot(jdir);
        
        minlnb_totSAVE=minlnb_tot(jdir);
        maxlnb_totSAVE=maxlnb_tot(jdir);
        
        lnbbins_pos_totSAVE=lnbbins_pos_tot(jdir);
        lnbbins_neg_totSAVE=lnbbins_neg_tot(jdir);
        
        
        bins_vapSAVE=bins_vap;
        bins_totSAVE=bins_tot;
                        
        exname=[exdir 'LNB_bel5ppmv.mat'];
        
        save(exname,'bins_vapSAVE','bins_totSAVE','me_lnb_abv_vapSAVE','me_lnb_bel_vapSAVE','minlnb_vapSAVE','maxlnb_vapSAVE','lnbbins_pos_vapSAVE','lnbbins_neg_vapSAVE',...
            'me_lnb_abv_totSAVE','me_lnb_bel_totSAVE','minlnb_totSAVE','maxlnb_totSAVE','lnbbins_pos_totSAVE','lnbbins_neg_totSAVE');
         

    end
    
    if saveflag(22)==1
		icediagsRADSAVE=icediagsRAD(jdir);
		exname=[exdir 'diagsRAD.mat'];
		save(exname,'icediagsRADSAVE');
    end
    
     if saveflag(23)==1
		distSAVE=dist(jdir);
		exname=[exdir 'fall_dist.mat'];
		save(exname,'distSAVE');
    end
    
    if saveflag(24)==1
		icediagstempSAVE=icediagsTEMP(jdir);
		exname=[exdir 'diagsTEMP.mat'];
		save(exname,'icediagstempSAVE');
    end
    
    if saveflag(25)==1
        
%         rhopertTimHSAVE=rhopertTimH;
%         rhopertTimHmaxSAVE=rhopertTimHmax;
%         rhopertTimHminSAVE=rhopertTimHmin;
%         
%         rhopertTimH2SAVE=rhopertTimH2;
%         rhopertTimHmax2SAVE=rhopertTimHmax2;
%         rhopertTimHmin2SAVE=rhopertTimHmin2; %mean temp pert
%         
%         rho5ppmv_vapSAVE=rho5ppmv_vap;
%         rho5ppmv_totSAVE=rho5ppmv_tot;
%         
%         rho5ppmv_vapnegSAVE=rho5ppmv_vapneg;
%         rho5ppmv_vapposSAVE=rho5ppmv_vappos;

		exname=[exdir 'rhopert5ppmv.mat'];
        savecom='save(exname';
        varlist={'rhopertTimH','rhopertTimHmax','rhopertTimHmin',...
                'rhopertTimH2','rhopertTimHmax2','rhopertTimHmin2',...
                'rho5ppmv_vap','rho5ppmv_tot',...
                'rho5ppmv_vapneg','rho5ppmv_vappos',...
                'rho5ppmv_totneg','rho5ppmv_totpos',...
                'low_prctiles','mid_prctiles','upp_prctiles','ztr_prctiles',...
                'tra','potemp_vap','potemp_tot'};
        
         varlist={'rhopertTimH','rhopertTimHmax','rhopertTimHmin',...
                'rhopertTimH2','rhopertTimHmax2','rhopertTimHmin2',...
                'rho5ppmv_vap','rho5ppmv_tot',...
                'rho5ppmv_vapneg','rho5ppmv_vappos',...
                'rho5ppmv_totneg','rho5ppmv_totpos',...
                'low_prctiles',...
                'tra','potemp_vap','potemp_tot'};
        
        for i=1:length(varlist)
            savestr{i}=[varlist{i} 'SAVE'];
            comm=[ varlist{i} 'SAVE=' varlist{i} '(jdir)'];
            eval(comm);     
			savecom=[savecom ',''' savestr{i} '''']        
        end

    savecom=[savecom ');'];
    
    eval(savecom);
        
%         savestr=savevar(savestr,'rho5ppmv_totneg');
%         savestr=savevar(savestr,'rho5ppmv_totpos');
%         
%         rho5ppmv_totneg(1).r(ikm,jj)=mean(rhopert(ikm,inon));
%         rho5ppmv_totpos(1).r(ikm,jj)=mean(rhopert(ikm,inon2));
% 
%         low_prctiles(1).t(1:size(TwoD.Q(:,:,10),1),jj,1:length(prcs))=(prctile(TwoD.Q(:,:,10)',prcs))';
%         mid_prctiles(1).t(1:size(TwoD.Q(:,:,11),1),jj,1:length(prcs))=(prctile(TwoD.Q(:,:,11)',prcs))';
%         upp_prctiles(1).t(1:size(TwoD.Q(:,:,12),1),jj,1:length(prcs))=(prctile(TwoD.Q(:,:,12)',prcs))';
%         ztr_prctiles(1).t(1:size(TwoD.Q(:,:,13),1),jj,1:length(prcs))=(prctile(TwoD.Q(:,:,13)',prcs))';
            
            
%             
% 		save(exname,'icediagstempSAVE');
    end
    
    if saveflag(26)==1
		exname=[exdir 'vaptotDists.mat'];
        varlist={'totdist','vapdist','meanvap','meantot'};
%        varlist={'vapdist'};
        savingCommands %saves variables in varlist (with SAVE on end)

    end
    
    if saveflag(262)==1
		exname=[exdir 'vaptotDists.mat'];
%        varlist={'totdist','vapdist','meanvap','meantot'};
        varlist={'vapdist','meanvap'};
        savingCommands %saves variables in varlist (with SAVE on end)

    end
    
    if saveflag(27)==1
        exname=[exdir 'windW_1-20_short.mat'];
        varlist={'TwoD_alltim'};
        
        savingCommands %saves variables in varlist (with SAVE on end)
    end
    
    if saveflag(28)==1
        exname=[exdir 'signTpert_lnb'];
        varlist={'tpertsign_lnb'};
        
        savingCommands %saves variables in varlist (with SAVE on end)
    end
    
    if saveflag(29)==1
        exname=[exdir 'qprctiles.mat'];
        varlist={'q_prctiles'};
        
        savingCommands %saves variables in varlist (with SAVE on end)
    end
    
    if saveflag(30)==1
        exname=[exdir 'dump62icediags'];
        varlist={'dump62icediags'};
        
        savingCommands %saves variables in varlist (with SAVE on end)
    end
    
    if saveflag(31)==1
        exname=[exdir 'TRACERacc3'];
        varlist={'TRACERacc3','CONDacc3','RAINacc3','LWCacc3','rhoacc3'};
        savingCommands %saves variables in varlist (with SAVE on end)
    end
    
    if saveflag(32)==1
        exname=[exdir 'INCgt1e8'];
        varlist={'rhoacc_alltim','INCacc_alltim','INCmaxacc_alltim','IWCacc_alltim'};
        savingCommands %saves variables in varlist (with SAVE on end)
    end
    
     if saveflag(33)==1
        exname=['c:/documents and settings/login/my documents/hibiscus/baloonDATA/balloon_diags'];
        varlist={'dmi','sdla','saw'};
        savingCommands %saves variables in varlist (with SAVE on end)
    end

    if saveflag(34)==1
         exname=[exdir 'rho_mean.mat'];
        varlist={'rho_prof'};
        savingCommands %saves variables in varlist (with SAVE on end)
    end
    
    if saveflag(35)==1
         exname=[exdir 'radar_10dbz.mat'];
        varlist={'n10dbz','n10dbz2'};
        savingCommands %saves variables in varlist (with SAVE on end)
    end
    
    if saveflag(36)==1
         exname=[exdir 'radar_echotops.mat'];
        varlist={'ntop'};
        savingCommands %saves variables in varlist (with SAVE on end)
    end
    
    if saveflag(37)==1
         exname=[exdir 'Tperts_all.mat'];
        varlist={'tpertTimH_full'};
        savingCommands %saves variables in varlist (with SAVE on end)
    end
    
    if saveflag(38)==1
         exname=[exdir 'CAPE_etc.mat'];
        varlist={'cape_etc'};
        savingCommands %saves variables in varlist (with SAVE on end)
    end
    
    if saveflag(39)==1
         exname=[exdir 'turbq_flux.mat'];
        varlist={'vappos_up',
                'vapneg_up',
                'vappos_d',
                'vapneg_d',
                'totpos_up',
                'totneg_up',
                'totpos_d',  %totpod_d and vappos_d might be missing
                'totneg_d' };
                    
        savingCommands %saves variables in varlist (with SAVE on end)
    end
        
   if saveflag(40)==1
        exname=[exdir 'vapfull.mat'];
        varlist={'vapfull'};
        savingCommands %saves variables in varlist (with SAVE on end)
    end 
    
    if saveflag(42)==1
        exname=[exdir 'fields_for_supersat_hslices2.mat'];
        varlist={'pressure_hslice','vap_hslice','potemp_hslice'};
        savingCommands %saves variables in varlist (with SAVE on end)
    end 
    
    if saveflag(43)==1
        exname=[exdir 'ice_hslices_justice.mat'];
        varlist={'iceMR_hslice','iceNC_hslice'}; %,'snowMR_hslice','snowNC_hslice'};
        savingCommands %saves variables in varlist (with SAVE on end)
    end 
    
    if saveflag(44)==1
        exname=[exdir 'lwc_width.mat'];
        varlist={'lwc_width'}; %,'snowMR_hslice','snowNC_hslice'};
        savingCommands %saves variables in varlist (with SAVE on end)
    end 
    
    if saveflag(45)==1
        exname=[exdir 'totice_44.mat'];
        varlist={'totice_44'}; %,'snowMR_hslice','snowNC_hslice'};
        savingCommands %saves variables in varlist (with SAVE on end)
    end 
    
    if saveflag(46)==1
        exname=[exdir 'dT.mat'];
        varlist={'dT_conv','dT_nonconv','dT_bubble'}; 
        savingCommands %saves variables in varlist (with SAVE on end)
    end 
    if saveflag(47)==1
        exname=[exdir 'tperts2.mat'];
        varlist={'tpertTimH_full','tpert_max','vappert_max','vappertTimH_full','capes'}; 
        savingCommands %saves variables in varlist (with SAVE on end)
    end     
        
    
    
end

fprintf(1,'\nSave program completed');