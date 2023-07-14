clear direc topdir

comp='lacieLap';
comp='uni';
%comp='docsLap';

ndir=120;

machineSTORE=[1 1 1 1 1 1 1 2 2 2 2 1 1 3 3 3 3 3 1 1 1 1 3 3 1 6 3 2 3 3 2 3 3 1 ...
        1 1 3 1 3 1 3 3 3 3 99 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3]; 
machineSTORE(732)=machineSTORE(73);

%execute these as default to set machinSTORE values
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
        topdir(12).dir='g:/runs/';
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
        topdir(26).dir='z:/';
        topdir(30).dir='g:/runs/';
        topdir(31).dir='g:/runs/';
        topdir(32).dir='g:/runs/';
        topdir(33).dir='g:/runs/';
        topdir(34).dir='g:/runs/';
        topdir(35).dir='g:/runs/';
        topdir(36).dir='g:/runs/';
        topdir(37).dir='g:/runs/';
        topdir(38).dir='g:/runs/';
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
        topdir(73).dir='g:/runs/';    %'z:/';              
        topdir(74).dir='g:/runs/';

        topdir(76).dir='z:/'; machineSTORE(length(topdir))=3;
        topdir(77).dir='z:/'; machineSTORE(length(topdir))=3; %3=Horace
        topdir(78).dir='g:/runs/'; machineSTORE(length(topdir))=3;
        topdir(79).dir='z:/'; machineSTORE(length(topdir))=3; %3=Horace
        topdir(80).dir='z:/'; machineSTORE(length(topdir))=3; %3=Horace
		topdir(81).dir='z:/'; machineSTORE(length(topdir))=3; %3=Horace
        topdir(82).dir='z:/'; machineSTORE(length(topdir))=3; %3=Horace
        topdir(83).dir='z:/'; machineSTORE(length(topdir))=3; %3=Horace    
        topdir(84).dir='z:/'; machineSTORE(length(topdir))=3; %3=Horace  
        topdir(85).dir='z:/'; machineSTORE(length(topdir))=3; %3=Horace 
        topdir(86).dir='z:/'; machineSTORE(length(topdir))=3; %3=Horace 
        topdir(87).dir='z:/'; machineSTORE(length(topdir))=3; %3=Horace 
        topdir(88).dir='z:/'; machineSTORE(length(topdir))=3; %3=Horace 
        topdir(89).dir='z:/'; machineSTORE(length(topdir))=3; %3=Horace 
        topdir(90).dir='z:/'; machineSTORE(length(topdir))=3; %3=Horace 
        topdir(91).dir='z:/'; machineSTORE(length(topdir))=3; %3=Horace 
        topdir(92).dir='z:/'; machineSTORE(length(topdir))=3; %3=Horace 
        topdir(93).dir='z:/'; machineSTORE(length(topdir))=3; %3=Horace 
        topdir(94).dir='z:/'; machineSTORE(length(topdir))=3; %3=Horace 
        
        for irun=95:101;
            topdir(irun).dir='z:/NM_runs/'; machineSTORE(irun)=3;
        end
        for irun=102:120;
            topdir(irun).dir='z:/'; machineSTORE(irun)=3;
        end
        for irun=105:108;
            topdir(irun).dir='z:/NM_runs/2e-3_ccn/'; machineSTORE(irun)=3;
        end
        
        for irun=[109:200];
            topdir(irun).dir='z:/'; machineSTORE(irun)=3;
        end
        
        for irun=[73 112 113 118 129 130 137 141:200]; %have put 129 and 130 onto y:
            topdir(irun).dir='y:/LEM/'; machineSTORE(irun)=3;
        end
        
        for irun=[133]; %have put 129 and 130 onto y:
            topdir(irun).dir='g:/runs_vapour_paper/'; machineSTORE(irun)=3;
        end
        
        for irun=[1:171]; %have put 129 and 130 onto y:
            add_ground_heights_store(irun).h=620/1000; %beeds to be km
        end
        
        for irun=[172:300]; %have put 129 and 130 onto y:
            topdir(irun).dir='y:/mmocca/'; machineSTORE(irun)=3;
            add_ground_heights_store(irun).h=0;
        end
                
        topdir(732).dir='z:/';  
        
switch comp     
case 'lacieLap'
    for i=1:ndir
        topdir(i).dir='e:/les/';
    end
    
case 'docsLap'
    for i=1:ndir
        topdir(i).dir='c:/documents and settings/G/my documents/lem/';
    end    
    
end
        
run_name(1).nam='1km res';
run_name(2).nam='1km res CCN 960 cm^{-3}';
run_name(3).nam='c:/cygwin/home/user/runs/';
run_name(4).nam='g:/runs/';
run_name(5).nam='h:/';
run_name(6).nam='g:/runs/';
run_name(7).nam='h:/';
run_name(8).nam='g:/runs/';
run_name(9).nam='';
run_name(10).nam='';
run_name(11).nam='';
run_name(12).nam='g:/runs/';
run_name(13).nam='h:/';
run_name(14).nam='250m res';
run_name(15).nam='3d';
run_name(16).nam='500m res';
run_name(17).nam='1km res';
run_name(18).nam='250m res, low updraught';
run_name(19).nam='c:/cygwin/home/user/runs/';
run_name(20).nam='c:/cygwin/home/user/runs/';
run_name(21).nam='c:/cygwin/home/user/runs/';
run_name(22).nam='c:/cygwin/home/user/runs/';
run_name(23).nam='250m res CCN 960cm^{-3}';
run_name(24).nam='c:/cygwin/home/user/runs/';
run_name(25).nam='Miles City 26.8';
run_name(26).nam='EMERALD Hector (Paul)';
run_name(27).nam='224m vert res';
run_name(28).nam='1km res ccn 960';   %done on HPCx and so wrong dynamics - updraughts too high
run_name(29).nam='1km res no rad';
run_name(30).nam='ccn 960 2000 km';
run_name(31).nam='HPCx_1kmccn960';
run_name(32).nam='10secs 8km damp';
run_name(33).nam='correct ice acc';
run_name(34).nam='Miles City 13.4 TH';        
run_name(35).nam='Miles City 13.4 TH Narrow';        
run_name(36).nam='Miles City 13.4 TH Narrow, dry lower';        
run_name(37).nam='2000 km';              
run_name(38).nam='13.02_2_13.4th';        
run_name(39).nam='4000km';        
run_name(40).nam='Miles City CCN=960 cm^{-3}';        
run_name(41).nam='2000km 20 mins';    run_name(41).nam='1 km resolution';    
run_name(42).nam='hightop';        
run_name(43).nam='2000km ccn960';        
run_name(44).nam='2000km xr=3.5km';        
run_name(45).nam= '20051116';       
run_name(46).nam='2000km 13.4th';       
run_name(47).nam='z_150'; 
run_name(48).nam='2km 2000km'; 
run_name(49).nam='2km 4000km';
run_name(50).nam='4000km_2';
run_name(51).nam='3d';
run_name(52).nam='2000km 6.7th';
run_name(53).nam='2000km 6.7th_noqv';
run_name(54).nam='1.25km';
run_name(55).nam='Redel_20mins';
run_name(56).nam='6.7th_1.67qv';
run_name(57).nam='highRH';
run_name(58).nam='nodry';
run_name(59).nam='3d2km_30red';
run_name(60).nam='tgraup2';
run_name(61).nam='4000km_3.5km';
run_name(62).nam='convap_allmoist';
run_name(63).nam='6.7th_20mins_3.5km';
run_name(64).nam='nodry_highRH';
run_name(65).nam='moist_noqv';
run_name(66).nam='3d2km_30';
run_name(67).nam='highT_test';
run_name(68).nam='highT_1.5';
run_name(69).nam='highT1.5_6.7th';
run_name(70).nam='noshear';
run_name(71).nam='2km_20mins';  run_name(71).nam='2D';
run_name(72).nam='6.7th_20mins';
run_name(73).nam='3D'; %300km case
run_name(732).nam='3D';
run_name(74).nam='CCN=960 cm^{-3}';

run_name(76).nam='3d 6.7th 1.67qv';
run_name(77).nam='3d ccn960';
run_name(78).nam='3D-med';  %'3d HR=6.7';
run_name(79).nam='3d 75km';
run_name(80).nam='NM test noheat';
run_name(81).nam='NM_hotspot_7km';
run_name(82).nam='NM_hotspot_1440LT';
run_name(83).nam='NM_moist_hotspot';
run_name(84).nam='NM_1440_nohot';
run_name(85).nam='NM_all_moist';
run_name(86).nam='NM_all_moist';
run_name(87).nam='NM_low_flux';

run_name(88).nam='150km_2.23th';
run_name(89).nam='150km_2.5qv';

run_name(90).nam='Montana 1e-3';
run_name(91).nam='NM_1100_nohot';
run_name(92).nam='NM_lowf_1e-3_15km';
run_name(93).nam='NM_lowf_1e-3_20km';

run_name(94).nam='new_coldpool';

run_name(95).nam= 'NM_ccn240/'; 
run_name(96).nam= 'NM_ccn480/'; 
run_name(97).nam= 'NM_ccn720/'; 
run_name(98).nam= 'NM_ccn960/'; 
run_name(99).nam= 'NM_newCCN_test/';

run_name(100).nam= 'Montana 2e-3/'; 
run_name(101).nam= 'Montana 4e-3/'; 

run_name(102).nam='newcp_6.7th';
run_name(103).nam='newcp_6.7th_10s';
run_name(104).nam='newcp_6.7th_1hr';

run_name(105).nam= 'CCN=240 cm^{-3}'; 
run_name(106).nam= 'CCN=480 cm^{-3}'; 
run_name(107).nam= 'CCN=720 cm^{-3}'; 
run_name(108).nam= 'CCN=960 cm^{-3}'; 

run_name(109).nam= '2000km_T1/'; 

run_name(110).nam= '2000km_26.8th/'; 
run_name(111).nam= '2000km_33.5th/'; 

run_name(112).nam= '3D-weak';    %'3d HR=2.23'; 
run_name(113).nam= '3d HR=4.46'; 

run_name(114).nam= '2000km_1e-3_1hour/';
run_name(115).nam= '2000km_1e-3_2hour/';
run_name(116).nam= '2000km_1e-3_3hour/';

run_name(117).nam= 'toga test ';
run_name(118).nam= '3D 150km';    %'150km trip20.1';

run_name(119).nam= '2000km 1e-3 0.5hour';
run_name(120).nam= 'toga test 304';
run_name(121).nam= 'toga test 306';

run_name(122).nam= 'toga test 1deg ';
run_name(123).nam= 'toga ad ';
run_name(124).nam= 'texas';
run_name(125).nam= 'new mexico';

run_name(126).nam= '150km ccn960 20.1th';

run_name(127).nam= 'texas superad';
run_name(128).nam= 'texas superad 1e-3';
run_name(129).nam= 'newmexico superad';
run_name(130).nam= 'newmexico superad 1e-3';
run_name(131).nam= 'toga 1e-3 1hr';
run_name(132).nam= 'toga 1deg 1e-3';

%don't put underscores in these!!
run_name(133).nam= 'CCN=960 cm^{-3}';    %'150km ccn960 non trip';
run_name(134).nam= '150km 20.1th';

run_name(135).nam= 'rhog 261';
run_name(136).nam= 'rhog 550';

run_name(137).nam= '75km bbigg50';
run_name(138).nam= '75km qv1.67';
run_name(139).nam= '75km qv3.34';
run_name(140).nam= '75km 2.23 ccn';

run_name(141).nam= '75km bbigg=0';
run_name(142).nam= '75km orig sounding';

run_name(143).nam= '75km IN 0.1';
run_name(144).nam= '75km HM 0.1';
run_name(145).nam= '75km HM+IN 0.1';

run_name(146).nam= '75km orig sounding 10mins';

run_name(147).nam= '75km 1km res';
run_name(148).nam= '75km orig hotpool';

run_name(149).nam= 'New Mexico supress';
run_name(150).nam= 'New Mexico supress 0.25e-3th';
run_name(151).nam= 'New Mexico supress 0.05e-3th';
run_name(152).nam= 'New Mexico supress 0.25e-3th 3.5km';
run_name(153).nam= 'New Mexico supress noth';
run_name(154).nam= 'New Mexico supress norh';

run_name(155).nam= 'New Mexico 2 deg supress norh';

run_name(156).nam= 'New Mexico 2 deg supress hightop';
run_name(157).nam= 'New Mexico 4 deg supress hightop';

run_name(158).nam= 'New Mexico 2 deg 10 mins';
run_name(159).nam= 'New Mexico 4 deg 40 mins';

run_name(160).nam= '75 km, 1km res';
run_name(161).nam= '75 km, orig';
run_name(162).nam= '75 km, 1km res 1.0th';
run_name(163).nam= '75 km, 1km res 2.23th'; 

run_name(164).nam= 'nm suppress 2deg ccn60';
run_name(165).nam= 'nm suppress 2deg ccn120';
run_name(166).nam= 'nm suppress 2deg ccn10';

run_name(167).nam= '150km ccn960 homog only';

run_name(168).nam= '75km 1km res 4.46th';
run_name(169).nam= '75km 1km res 8.92th';

run_name(170).nam= 'NM suppress 3d';

run_name(171).nam= 'NM suppress 2deg ccn960';

run_name(172).nam= 'Texas supp 4deg 1e-3th';
run_name(173).nam= 'Texas supp 4deg 2e-3th';
run_name(174).nam= 'Texas supp 4deg 1e-3th 14km';
run_name(175).nam= 'Texas supp 4deg 1e-3th 21km';
run_name(176).nam= 'Texas supp 4deg 1e-3th bigger dom';
run_name(177).nam= 'Texas supp 4deg 0.5e-3th';
run_name(178).nam= 'Miles supp 4deg 1e-3th';
run_name(179).nam= 'Texas supp 4deg 1e-3th Hector';
run_name(180).nam= 'Texas supp4deg 1e-3 Hector_fastsse';

run_name(181).nam= 'Texas supp4deg 0.25e-3th';
run_name(182).nam= 'Texas supp4deg 0.1e-3th';
run_name(183).nam= 'Miles supp4deg 0.05e-3th';
run_name(184).nam= 'Miles supp4deg 0.1e-3 Hector_fastsse';
run_name(185).nam= 'Miles supp4deg 0.2e-3 Hector_fastsse';

run_name(186).nam= 'Texas 10mins 0.1e-3';
run_name(187).nam= 'Texas supp4deg 1e-5th';
run_name(188).nam= 'Texas supp4deg noheat';
run_name(189).nam= 'Texas supp6deg 20mins';

run_name(190).nam= 'tx_20mins_6deg';
run_name(191).nam= 'Texas 1e-5 rh0';
run_name(192).nam= 'Miles 4deg 0.4 rh';
run_name(193).nam= 'Miles 6deg';
run_name(194).nam= 'Miles4_tau60';
run_name(195).nam= 'Texas_noheat_irh0';
run_name(196).nam= 'Texas10min_supp6';
run_name(197).nam= 'tx_1e-5_irh0';
run_name(198).nam= 'Texas supp6 tau60';
run_name(199).nam= 'Texas supp6 tau60 3.5km';
run_name(200).nam= 'Miles supp4 tau60 3.5km';
run_name(201).nam= 'Texas supp6 tau60 3.5km 10mins';
run_name(202).nam= 'Texas supp6 tau60 3.5km 10mins 1e-5';

run_name(203).nam= 'Miles 6deg 120 CCN';
run_name(204).nam= 'Miles 6deg 360 CCN';
run_name(205).nam= 'Miles 6deg 480 CCN';
run_name(206).nam= 'Miles 6deg 600 CCN';
run_name(207).nam= 'Miles 6deg 720 CCN';
run_name(208).nam= 'Miles 6deg 840 CCN';
run_name(209).nam= 'Miles 6deg 960 CCN';

run_name(210).nam= 'Texas supp 8';

run_name(211).nam= 'Miles 6deg Vgraup';
run_name(212).nam= 'Miles 6deg MidLat graup';

run_name(213).nam= 'MilesV 6deg 120 CCN';
run_name(214).nam= 'MilesV 6deg 360 CCN';
run_name(215).nam= 'MilesV 6deg 480 CCN';
run_name(216).nam= 'MilesV 6deg 600 CCN';
run_name(217).nam= 'MilesV 6deg 720 CCN';
run_name(218).nam= 'MilesV 6deg 840 CCN';
run_name(219).nam= 'MilesV 6deg 960 CCN';


run_name(220).nam= 'tx 6deg 120 CCN';
run_name(221).nam= 'tx 6deg 360 CCN';
run_name(222).nam= 'tx 6deg 480 CCN';
run_name(223).nam= 'tx 6deg 600 CCN';
run_name(224).nam= 'tx 6deg 720 CCN';
run_name(225).nam= 'tx 6deg 840 CCN';
run_name(226).nam= 'tx 6deg 960 CCN';
run_name(227).nam= 'tx 6deg 60 CCN';


run_name(228).nam= 'tx 6deg 60 CCN';
run_name(229).nam= 'tx 6deg 120 CCN';
run_name(230).nam= 'tx 6deg 240 CCN';
run_name(231).nam= 'tx 6deg 360 CCN';
run_name(232).nam= 'tx 6deg 480 CCN';
run_name(233).nam= 'tx 6deg 600 CCN';
run_name(234).nam= 'tx 6deg 720 CCN';
run_name(235).nam= 'tx 6deg 840 CCN';
run_name(236).nam= 'tx 6deg 960 CCN';

run_name(237).nam= 'tx morr mphys';

run_name(238).nam= 'tx morr mphys';
run_name(239).nam= 'tx morr homog';
run_name(240).nam= 'tx morr psacw';
run_name(241).nam= 'tx morr homog t0.01';
run_name(242).nam= 'tx morr homog lwc5.1';
run_name(243).nam= 'tx morr homog justLWC';
run_name(244).nam= 'tx morr psacw 0.35';
run_name(245).nam= 'tx morr psacw 0.85';

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

loadflag=zeros([1 150]); 
%diags to load

%loadflag([3 5 6 7 8 11 12 17 29])=1; %usual suite as needed from the save script

%loadflag([11])=1;

%loadflag([3])=1; %maxw

%loadflag([5 7])=1;  %dqtot and dqvap

% endt{1}='58'; %direc(1).dir='h:/5ppmv_qnums2/';
% endtALL{1}='58';

%loadflag([8 11])=1; %11=icediag_nums

%loadflag([8])=1; %8=icediagsALL
%loadflag([9])=1; %8=icediag - old skool microphysical rates (before the number change thing for homog freezing etc.)

%loadflag([6])=1; %4=vap&tot prctiles

%loadflag(5)=1; 

%loadflag(17)=1; %for max ice supersat time height
%loadflag(18)=1; %rhopert and temppert
%loadflag(19)=1; %rhopert and temppert
%loadflag(20)=1; %mean postive and neg LNB
%loadflag(21)=1; %mean lnb pos & neg in low vap and tot plus binned distributions
%loadflag(22)=1; %radiation stats
%loadflag(23)=1; %fall speed stats
%loadflag(24)=1; %temp data xxx_TH
%loadflag(25)=1; %density perturbations
%loadflag(26)=1; %vapour and tot distributions
%loadflag(27)=1; % all updraught winds
%loadflag(29)=1; %q-prctiles
%loadflag(30)=1; %icediagsALL for just dump 62 (in order to get ACC_A for eq. model)
%loadflag(31)=1; %tracerACC diags

%loadflag(32)=1; %INCgt1e8
%loadflag(33)=1; %updraught width
%loadflag(34)=1; %SF4 soundings
%loadflag(35)=1; %max tracer time-height
%loadflag(36)=1; %rho mean
%loadflag(37)=1; %length of domain with > n dBz radar echo at each height n10dbz
%loadflag(38)=1; %length of domain with > n dBz radar echo at each height
%loadflag(390=1; %temperature perts of full domain  %these available for the high CCN case
%loadflag(40)=1; %vappos etc. w'qv' separated into updraughts and down and positive and neg qv from 2-d fields
%loadflag(41)=1; %full 2d vapour fields 
%loadflag(42)=1; %horizontal slice at h=17.9km of fields for supersaturation calc (temp,press,vap)
%loadflag(43)=1; %horizontal slice at h=17.9km of ice no and mass (just ice, not snow or graupel)
%loadflag(45)=1; %tot ice or water (dump 44) tot_ice44

%loadflag(46)=1; %dT_conv, etc
%loadflag(47)=1; %tpert_max etc.



%*******************************************************************************************************
justname=1; %flag to say just want the name for loading purposes (ignores loadflags)
gridflag=1; %flag to say whether to load in Grid
  
ishift=0; %flag to say that want icediagsALL and dq_tot, etc. data in 2000km cases shifted accross in time by 2 indices
%%%%%%%% make sure to set to zero if are loading in order to contiue processing ****


%loadselect=[14 16 17 18]; %files to load
%loadselect=[17];
%loadselect=[17 23 1 14];
%loadselect=[14 23 24];
loadselect=[14 16 17]; %diff res runs
%loadselect=[25 34 33 35 14]; 
%loadselect=[1]; %dmiDamp_2
loadselect=[1 2 14 23]; %1km, ,250m norm annd high ccn
loadselect=[14 23]; %1km, ,250m norm annd high ccn

% loadselect=[20 14 16:18]; %
% loadselect=[17 37 39];
% loadselect=[37 41 44 48 49 50 52];
% loadselect=[37 41 44 52 62 57 58 56 65];
% loadselect=[62 57 58 56 65];
% loadselect=[67 62 52 41];
% loadselect=[67 69 52 44 68];
% loadselect=[44 48 52 37 41];
% loadselect=[37 48 44 52 41];
% loadselect=[41 74]; %low/high CCN
loadselect=[73 71]; %3D 2D
%loadselect=[118 71]; %3D 2D

%loadsele=[78 79];
%loadselectct=[41 71]; %1km 2km
%loadselect=[34 35 36]; %Miles City cases
%loadselect=[79 73 78];
%loadselect=[732 79];
%loadselect=[81:84 80 85 86 87 90 91];
%loadselect=[84 90];
%loadselect=[78];
%loadselect=[732 78 88 89 71];


loadselect=[112];
%loadselect=[73];


%loadselect=[90 100 101];
%loadselect=[95:98];
%loadselect=[105:108];
%loadselect=[110 111 71 72 732];
%loadselect=[111 71 732 78];
%loadselect=[112 732 78];
%loadselect=[117 120:123 131 132];
%loadselect=[125 129 130];
%loadselect=[90 100 101];
%loadselect=[105:108];

%loadselect=[73 133]; % 3D_2km_20mins, 150km_ccn960_non_trip
%loadselect=[118 133]; % 150km_trip20.1, 150km_ccn960_non_trip

%loadselect=[732 118];

%loadselect=[732 78 112]; 
%loadselect=[71 73 78 113 112]; 
%loadselect=[73 78 112]; %control, medium and weak cases 

%loadselect=[156 170 157];  %emm comparison - 2deg suppress hightop, w/ diags
%loadselect=[156 171];  %emm comparisons

%loadselect=[156 164 165 166]; 

%loadselect=[133]; 

%loadselect=[112 113 118 71]; 

%loadselect=[118 112 113 133]; %strong, medium, weak and high CCN cases
%loadselect=[73 112 113 133]; %strong, medium, weak and high CCN cases

%loadselect=[118 133 167];  %3d and 3d homogen onl%stratoy run (run where ccn effects only apply to the ice numbers at homog

%loadselect=[162 163 168 169 118]; %1km res runs

%loadselect=[73 118 137]; %3D 20.1th runs with 300 and 150 km and 75 km domain sizes

%loadselect=[148 73]; %75km_orig_hotpool and 300km 3D case 
%loadselect=[118 169 73 160]; %75km_1km_8.92th

%loadselect=[193 211 212];
%loadselect=[201 202 210];

%loadselect=[176:177 181:182];
%loadselect=[203 193 204:209 210];

%loadselect=[213 211 214:219]; %ccn sens for miles city with vgraup
%[27,44,27,27,28,22,23,27]

%loadselect=[227 220 210 221:226]; %ccn sens for texas case
%loadselect=[243 244 245]; %
%loadselect=[230 239]; %normal and homog frist case (0.05s tstep and 1.7e-3 LWC)

loadselect=[238 240:245 230];
loadselect=[230 244];

%[33 33 33 33 33 33 33 33 33]



%**********************************************************************************************************

        %direc(1).dir='c:/documents and settings/tom/my documents/dan/dmi1715_5ppmv_25km_3/';
%direc(2).dir='c:/documents and settings/tom/my documents/dan/dmidamp_2/';
direc(1).dir=[topdir(1).dir 'dmidamp_2/']; %zdmp=22.858km
direc(2).dir=[topdir(2).dir 'damp_ccn960/']; %zdmp=22.858km  
%direc(4).dir='c:/documents and settings/tom/my documents/dan/dmi1715_5ppmv/';
direc(3).dir=[topdir(3).dir 'damp_inx10/']; %zdmp=22.858km
direc(4).dir=[topdir(4).dir '500mres/']; %zdmp=22.858km

direc(5).dir=[topdir(5).dir '5ppmv_qnums2/']; %xc for the hotpool =0 so think didn't run properly **zdmp=27.946km

%direc(1).dir='c:/cygwin/home/user/runs/5ppmv/';
direc(6).dir=[topdir(6).dir '5ppmv/']; %no icediag_nums   **zdmp=27.946km

direc(7).dir=[topdir(7).dir '5ppmv_1000km_3/']; %has icediag_nums  **zdmp=27.946km

direc(8).dir=[topdir(8).dir '5ppmv_ccn960_2/']; %Bezier runs  **zdmp=27.946km %500km run
direc(9).dir=[topdir(9).dir '5ppmv_indiv10/']; %  **zdmp=27.946km                  "
direc(10).dir=[topdir(10).dir '5ppmv_th26.8_2/']; % **zdmp=27.946km                "
direc(11).dir=[topdir(11).dir '4ppmv_th26.8_2/']; % **zdmp=27.946km           %500km run  
 
direc(12).dir=[topdir(12).dir '4ppmv/']; % zdmp=27.946km                      %500km run 

direc(13).dir=[topdir(13).dir 'lowE/']; % **zdmp=27.946km

direc(14).dir=[topdir(14).dir '250mres_1000km/']; % **zdmp=22.858km %no ozone tracer (14 and 15 defunct here and should be sorted or removed)
direc(15).dir=[topdir(15).dir '3dNewton/']; 
direc(16).dir=[topdir(16).dir '500mNewt/']; %1000km 
direc(17).dir=[topdir(17).dir '1kmNewt/']; %1000km
direc(18).dir=[topdir(18).dir '250m13.4Newt/']; %1000km

direc(19).dir=[topdir(19).dir '13.02.12utc/']; %13th Feb case study - started at 15UTC (not 12utc but 12 LT). 250km domain but warmpool put at 
                                                %250km so will only get half of it in domain as was using uncorrected code
direc(20).dir=[topdir(20).dir '13.02_2/']; %as above but with warmpool at x=10km so all of it will be within domain - both done on PC

direc(21).dir=[topdir(21).dir '13.02_inx10/']; 
direc(22).dir=[topdir(22).dir '13.02_ccn960/'];

direc(23).dir=[topdir(23).dir 'ccn960Newt/']; 
direc(24).dir=[topdir(24).dir 'inx10newt/'];

direc(25).dir=[topdir(25).dir 'MilesCity/'];

direc(26).dir=[topdir(26).dir]; %Paul's EMERALD CCN case given to Stewart for Hector sounding - need to change RUN001 to RUN003 for control
direc(26).dir=['H:/hector_Paul/']; %where diags are saved - runs stored in above
direc(27).dir=[topdir(25).dir 'z_150/']; %where diags are saved - runs stored in above
direc(28).dir=[topdir(28).dir '1kmccn960_HPCx/']; %from HPCx, 4 processors
direc(29).dir=[topdir(29).dir '1km_norad/']; %from Horace, 2 processors
direc(30).dir=[topdir(30).dir 'ccn960_2000km/']; %from Horace (machine=3), 2 processors  %no ozone tracer (14 and 15 defunct here and should be sorted or removed)
direc(31).dir=[topdir(31).dir 'HPCx_1kmccn960/']; %from HPCx, 4 processors
direc(32).dir=[topdir(32).dir 'damp_10secs_8km/']; 
direc(33).dir=[topdir(33).dir 'iceacc_1000km/']; 

direc(34).dir=[topdir(34).dir 'MilesCity13.4/'];
direc(35).dir=[topdir(35).dir 'Miles13.4_xr2000/'];
direc(36).dir=[topdir(36).dir 'Miles13.4_dry_bot/'];

direc(37).dir=[topdir(37).dir '2000km/']; %Horace run, 4 processors

direc(38).dir=[topdir(38).dir '13.02_2_13.4th/'];
direc(39).dir=[topdir(39).dir '4000km/']; %Horace run, 4 processors
direc(40).dir=[topdir(40).dir 'milescity13.4_ccn960/']; %Horace run, 4 processors

direc(41).dir=[topdir(41).dir '2000km_20mins/']; %Horace run, 2 processors
direc(42).dir=[topdir(42).dir 'hightop/']; %Horace run, 2 processors
direc(43).dir=[topdir(43).dir '2000km_ccn960_2/']; %Horace run, 2 processors
direc(44).dir=[topdir(44).dir '2000km_xr3500/']; %Horace run, 2 processors
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
direc(63).dir=[topdir(63).dir '6.7th_20mins_3.5km/'];
direc(64).dir=[topdir(64).dir 'nodry_highRH/'];
direc(65).dir=[topdir(65).dir 'moist_noqv/'];
direc(66).dir=[topdir(66).dir '3d2km_30/'];
direc(67).dir=[topdir(67).dir 'highT_test/'];
direc(68).dir=[topdir(68).dir 'highT_1.5/'];
direc(69).dir=[topdir(69).dir 'highT1.5_6.7th/'];
direc(70).dir=[topdir(70).dir 'noshear/'];
direc(71).dir=[topdir(71).dir '2km_20mins/'];
direc(72).dir=[topdir(72).dir '6.7th_20mins/'];
direc(73).dir=[topdir(73).dir '3d2km_30_20mins/'];   %''];
direc(732).dir=[topdir(732).dir ''];
direc(74).dir=[topdir(74).dir 'ccn960_20mins/'];

direc(76).dir=[topdir(76).dir '3d2km_30red/']; %6.7th and 1.67qv
direc(77).dir=[topdir(77).dir '3d_ccn960/']; 
direc(78).dir=[topdir(78).dir '3d2km_6.7th_20mins/']; 
direc(79).dir=[topdir(length(direc)).dir '3d_75km/']; 
direc(80).dir=[topdir(length(direc)).dir 'NM_tests_noheat/']; 
direc(81).dir=[topdir(length(direc)).dir 'NM_hotspot_7km/'];
direc(82).dir=[topdir(length(direc)).dir 'NM_hotspot_1440LT/']; 
direc(83).dir=[topdir(length(direc)).dir 'NM_moist_hotspot/']; 
direc(84).dir=[topdir(length(direc)).dir 'NM_1440_nohot/']; 
direc(85).dir=[topdir(length(direc)).dir 'NM_all_moist/']; 
direc(86).dir=[topdir(length(direc)).dir 'NM_no_inversion/']; 
direc(87).dir=[topdir(length(direc)).dir 'NM_low_flux/']; 

direc(88).dir=[topdir(length(direc)).dir '150km_2.23th/'];
direc(89).dir=[topdir(length(direc)).dir '150km_2.5qv/'];

direc(90).dir=[topdir(length(direc)).dir 'NM_lowf_1e-3/'];
direc(91).dir=[topdir(length(direc)).dir 'NM_1100_nohot/'];
direc(92).dir=[topdir(length(direc)).dir 'NM_lowf_1e-3_15km/'];
direc(93).dir=[topdir(length(direc)).dir 'NM_lowf_1e-3_20km/'];
direc(94).dir=[topdir(length(direc)).dir 'new_coldpool/'];

irun=95;  direc(irun).dir=[topdir(irun).dir 'NM_ccn240/']; 
irun=96;  direc(irun).dir=[topdir(irun).dir 'NM_ccn480/']; 
irun=97;  direc(irun).dir=[topdir(irun).dir 'NM_ccn720/']; 
irun=98;  direc(irun).dir=[topdir(irun).dir 'NM_ccn960/']; 
irun=99; direc(irun).dir=[topdir(irun).dir 'NM_newCCN_test/'];
irun=100;  direc(irun).dir=[topdir(irun).dir 'NM_lowf_2e-3/']; 
irun=101; direc(irun).dir=[topdir(irun).dir 'NM_lowf_4e-3/']; 



irun=102; direc(irun).dir=[topdir(irun).dir 'newcp_6.7th/']; 
irun=103; direc(irun).dir=[topdir(irun).dir 'newcp_6.7th_10s/']; 
irun=104; direc(irun).dir=[topdir(irun).dir 'newcp_6.7th_1hr/']; 

irun=105;  direc(irun).dir=[topdir(irun).dir 'NM_2e-3_240ccn/']; 
irun=106;  direc(irun).dir=[topdir(irun).dir 'NM_2e-3_480ccn/']; 
irun=107;  direc(irun).dir=[topdir(irun).dir 'NM_2e-3_720ccn/']; 
irun=108;  direc(irun).dir=[topdir(irun).dir 'NM_2e-3_960ccn/']; 
irun=109;  direc(irun).dir=[topdir(irun).dir '2000km_T1/']; 
irun=110;  direc(irun).dir=[topdir(irun).dir '2000km_26.8th/']; 
irun=111;  direc(irun).dir=[topdir(irun).dir '2000km_33.5th/']; 

irun=112;  direc(irun).dir=[topdir(irun).dir '150km_2.23th/']; 
irun=113;  direc(irun).dir=[topdir(irun).dir '150km_4.46th/']; 

irun=114;  direc(irun).dir=[topdir(irun).dir '2000km_1e-3_1hour/']; 
irun=115;  direc(irun).dir=[topdir(irun).dir '2000km_1e-3_2hour/']; 
irun=116;  direc(irun).dir=[topdir(irun).dir '2000km_1e-3_3hour/']; 

irun=117;  direc(irun).dir=[topdir(irun).dir 'toga_test/']; 
irun=118;  direc(irun).dir=[topdir(irun).dir '150km_trip20.1/']; 

irun=119;  direc(irun).dir=[topdir(irun).dir '2000km_1e-3_0.5hour/']; 

irun=120;  direc(irun).dir=[topdir(irun).dir 'toga_test_304/']; 
irun=121;  direc(irun).dir=[topdir(irun).dir 'toga_test_306/']; 

irun=122;  direc(irun).dir=[topdir(irun).dir 'toga_test_1deg/'];
irun=123;  direc(irun).dir=[topdir(irun).dir 'toga_ad/'];
irun=124;  direc(irun).dir=[topdir(irun).dir 'texas/'];
irun=125;  direc(irun).dir=[topdir(irun).dir 'new_mexico/'];

irun=126;  direc(irun).dir=[topdir(irun).dir '150km_ccn960_20.1th/'];

irun=127;  direc(irun).dir=[topdir(irun).dir 'texas_superad/'];
irun=128;  direc(irun).dir=[topdir(irun).dir 'texas_superad_1e-3/'];
irun=129;  direc(irun).dir=[topdir(irun).dir 'newmexico_superad/'];
irun=130;  direc(irun).dir=[topdir(irun).dir 'newmexico_superad_1e-3/'];
irun=131;  direc(irun).dir=[topdir(irun).dir 'toga_1e-3_1hr/'];
irun=132;  direc(irun).dir=[topdir(irun).dir 'toga_1deg_1e-3/'];

irun=133;  direc(irun).dir=[topdir(irun).dir '150km_ccn960_non_trip/'];
irun=134;  direc(irun).dir=[topdir(irun).dir '150km_20.1th/']; %actually is 6.7th - so can only do
%                       ccn comparison using 150km/300km cases but can't be sure of domain size effects
irun=135;  direc(irun).dir=[topdir(irun).dir 'rhog_261/']; %
irun=136;  direc(irun).dir=[topdir(irun).dir 'rhog_550/']; %

irun=137;  direc(irun).dir=[topdir(irun).dir '75km_bbigg50/']; %
irun=138;  direc(irun).dir=[topdir(irun).dir '75km_qv1.67/']; %
irun=139;  direc(irun).dir=[topdir(irun).dir '75km_qv3.34/']; %
irun=140;  direc(irun).dir=[topdir(irun).dir '75km_2.23_ccn/']; %
irun=141;  direc(irun).dir=[topdir(irun).dir '75km_bbigg0/']; %

irun=142;  direc(irun).dir=[topdir(irun).dir '75km_orig_sounding/']; %

irun=143;  direc(irun).dir=[topdir(irun).dir '75km_Meyers0.1/']; %
irun=144;  direc(irun).dir=[topdir(irun).dir '75km_HM_0.1/']; %
irun=145;  direc(irun).dir=[topdir(irun).dir '75km_Meyers0.1_HM_0.1/']; %

irun=146;  direc(irun).dir=[topdir(irun).dir '75km_orig_sounding_10mins/']; %

irun=147;  direc(irun).dir=[topdir(irun).dir '75km_1kmres/']; %
irun=148;  direc(irun).dir=[topdir(irun).dir '75km_orig_hotpool/']; %

irun=149;  direc(irun).dir=[topdir(irun).dir 'newmexico_supress/']; %
irun=150;  direc(irun).dir=[topdir(irun).dir 'nm_supress_0.25e-3th/']; %
irun=151;  direc(irun).dir=[topdir(irun).dir 'nm_suppress_5e-5th/']; %
irun=152;  direc(irun).dir=[topdir(irun).dir 'nm_suppress_3.5km/']; %
irun=153;  direc(irun).dir=[topdir(irun).dir 'nm_suppress_noth/']; %
irun=154;  direc(irun).dir=[topdir(irun).dir 'nm_suppress_norh/']; %
irun=155;  direc(irun).dir=[topdir(irun).dir 'nm_supp_2deg_norh/']; %
irun=156;  direc(irun).dir=[topdir(irun).dir 'nm_supp_2deg_hightop/']; %

irun=157;  direc(irun).dir=[topdir(irun).dir 'nm_supp_4deg/']; %

irun=158;  direc(irun).dir=[topdir(irun).dir 'nm_2deg_10mins/']; %
irun=159;  direc(irun).dir=[topdir(irun).dir 'nm_4deg_40mins/']; %

irun=160;  direc(irun).dir=[topdir(irun).dir '75km_1kmres2/']; %
irun=161;  direc(irun).dir=[topdir(irun).dir '75km_orig_18thJune07/']; %
irun=162;  direc(irun).dir=[topdir(irun).dir '75km_1kmres_1.0th/']; %
irun=163;  direc(irun).dir=[topdir(irun).dir '75km_1kmres_2.23th/']; %

irun=164;  direc(irun).dir=[topdir(irun).dir 'nm_supp_2deg_ccn60/']; %
irun=165;  direc(irun).dir=[topdir(irun).dir 'nm_supp_2deg_ccn120/']; %
irun=166;  direc(irun).dir=[topdir(irun).dir 'nm_supp_2deg_ccn10/']; %

irun=167;  direc(irun).dir=[topdir(irun).dir '150km_ccn960_homog_only/']; %

irun=168;  direc(irun).dir=[topdir(irun).dir '75km_1kmres_4.46th/']; %
irun=169;  direc(irun).dir=[topdir(irun).dir '75km_1kmres_8.92th/']; %
irun=170;  direc(irun).dir=[topdir(irun).dir 'nm3d_75km/']; %
irun=171;  direc(irun).dir=[topdir(irun).dir 'nmSep07_ccn960/']; %

irun=172;  direc(irun).dir=[topdir(irun).dir 'texas_supp4deg_1e-3/']; %
irun=173;  direc(irun).dir=[topdir(irun).dir 'texas_supp4deg_2e-3/']; %
irun=174;  direc(irun).dir=[topdir(irun).dir 'texas_supp4deg_xr14km/']; %
irun=175;  direc(irun).dir=[topdir(irun).dir 'texas_supp4deg_xr21km/']; %
irun=176;  direc(irun).dir=[topdir(irun).dir 'texas_supp4deg_1e-3_bigger_dom/']; %
irun=177;  direc(irun).dir=[topdir(irun).dir 'texas_supp4deg_0.5e-3/']; %
irun=178;  direc(irun).dir=[topdir(irun).dir 'miles_supp4deg_1e-3/']; %
irun=179;  direc(irun).dir=[topdir(irun).dir 'texas_supp4deg_1e-3_Hector/']; %
irun=180;  direc(irun).dir=[topdir(irun).dir 'texas_supp4deg_1e-3_Hector_fastsse/']; %

irun=181;  direc(irun).dir=[topdir(irun).dir 'texas_supp4deg_0.25e-3/']; %
irun=182;  direc(irun).dir=[topdir(irun).dir 'texas_supp4deg_0.1e-3/']; %
irun=183;  direc(irun).dir=[topdir(irun).dir 'miles_supp4deg_0.05e-3/']; %
irun=184;  direc(irun).dir=[topdir(irun).dir 'miles_supp4deg_0.1e-3/']; %
irun=185;  direc(irun).dir=[topdir(irun).dir 'miles_supp4deg_0.2e-3/']; %

irun=186;  direc(irun).dir=[topdir(irun).dir 'texas_10mins/']; %
irun=187;  direc(irun).dir=[topdir(irun).dir 'texas_1e-5/']; %
irun=188;  direc(irun).dir=[topdir(irun).dir 'texas_noheat/']; %
irun=189;  direc(irun).dir=[topdir(irun).dir 'tx_20mins_6deg/']; %

irun=190;  direc(irun).dir=[topdir(irun).dir 'tx_20mins_6deg2/']; %
irun=191;  direc(irun).dir=[topdir(irun).dir 'Tx_1e-5_irh0/']; %
irun=192;  direc(irun).dir=[topdir(irun).dir 'Miles_4deg_rh0.4/']; %
irun=193;  direc(irun).dir=[topdir(irun).dir 'Miles_6deg/']; %
irun=194;  direc(irun).dir=[topdir(irun).dir 'Miles4_tau60/']; %
irun=195;  direc(irun).dir=[topdir(irun).dir 'Texas_noheat_irh0/']; %
irun=196;  direc(irun).dir=[topdir(irun).dir 'Texas10min_supp6/']; %
irun=197;  direc(irun).dir=[topdir(irun).dir 'tx_1e-5_irh0/']; %
irun=198;  direc(irun).dir=[topdir(irun).dir 'tx_6deg_tau60/']; %
irun=199;  direc(irun).dir=[topdir(irun).dir 'tx_6deg_tau60_3.5km/']; %
irun=200;  direc(irun).dir=[topdir(irun).dir 'miles4_tau60_3.5km/']; %
irun=201;  direc(irun).dir=[topdir(irun).dir 'tx_6deg_tau60_3.5km_10mins/']; %
irun=202;  direc(irun).dir=[topdir(irun).dir 'tx_10mins_1e-5/'];

irun=203;  direc(irun).dir=[topdir(irun).dir 'miles6_ccn120/']; %
irun=204;  direc(irun).dir=[topdir(irun).dir 'miles6_ccn360/']; %
irun=205;  direc(irun).dir=[topdir(irun).dir 'miles6_ccn480/']; %
irun=206;  direc(irun).dir=[topdir(irun).dir 'miles6_ccn600/']; %
irun=207;  direc(irun).dir=[topdir(irun).dir 'miles6_ccn720/']; %
irun=208;  direc(irun).dir=[topdir(irun).dir 'miles6_ccn840/']; %
irun=209;  direc(irun).dir=[topdir(irun).dir 'miles6_ccn960/']; %

irun=210;  direc(irun).dir=[topdir(irun).dir 'tx_supp8/'];

irun=211;  direc(irun).dir=[topdir(irun).dir 'miles_vgraup/'];
irun=212;  direc(irun).dir=[topdir(irun).dir 'miles_midlat/'];

irun=213;  direc(irun).dir=[topdir(irun).dir 'milesv_ccn120/']; %
irun=214;  direc(irun).dir=[topdir(irun).dir 'milesv_ccn360/']; %
irun=215;  direc(irun).dir=[topdir(irun).dir 'milesv_ccn480/']; %
irun=216;  direc(irun).dir=[topdir(irun).dir 'milesv_ccn600/']; %
irun=217;  direc(irun).dir=[topdir(irun).dir 'milesv_ccn720/']; %
irun=218;  direc(irun).dir=[topdir(irun).dir 'milesv_ccn840/']; %
irun=219;  direc(irun).dir=[topdir(irun).dir 'milesv_ccn960/']; %


irun=220;  direc(irun).dir=[topdir(irun).dir 'tx8_ccn120/']; %
irun=221;  direc(irun).dir=[topdir(irun).dir 'tx8_ccn360/']; %
irun=222;  direc(irun).dir=[topdir(irun).dir 'tx8_ccn480/']; %
irun=223;  direc(irun).dir=[topdir(irun).dir 'tx8_ccn600/']; %
irun=224;  direc(irun).dir=[topdir(irun).dir 'tx8_ccn720/']; %
irun=225;  direc(irun).dir=[topdir(irun).dir 'tx8_ccn840/']; %
irun=226;  direc(irun).dir=[topdir(irun).dir 'tx8_ccn960/']; %
irun=227;  direc(irun).dir=[topdir(irun).dir 'tx8_ccn60/']; %

irun=228;  direc(irun).dir=[topdir(irun).dir 'tx8v_ccn60/']; %
irun=229;  direc(irun).dir=[topdir(irun).dir 'tx8v_ccn120/']; %
irun=230;  direc(irun).dir=[topdir(irun).dir 'tx8v_ccn240/']; %
irun=231;  direc(irun).dir=[topdir(irun).dir 'tx8v_ccn360/']; %
irun=232;  direc(irun).dir=[topdir(irun).dir 'tx8v_ccn480/']; %
irun=233;  direc(irun).dir=[topdir(irun).dir 'tx8v_ccn600/']; %
irun=234;  direc(irun).dir=[topdir(irun).dir 'tx8v_ccn720/']; %
irun=235;  direc(irun).dir=[topdir(irun).dir 'tx8v_ccn840/']; %
irun=236;  direc(irun).dir=[topdir(irun).dir 'tx8v_ccn960/']; %

irun=237;  direc(irun).dir=[topdir(irun).dir 'morr_test2/']; %
irun=238;  direc(irun).dir=[topdir(irun).dir 'morr_nomin_hector/']; %
irun=239;  direc(irun).dir=[topdir(irun).dir 'morr_homog/']; %
irun=240;  direc(irun).dir=[topdir(irun).dir 'morr_nomin_hector_2/']; %
irun=241;  direc(irun).dir=[topdir(irun).dir 'homog_t0.01/']; %
irun=242;  direc(irun).dir=[topdir(irun).dir 'homog_lwc5.1/']; %
irun=243;  direc(irun).dir=[topdir(irun).dir 'homog_justlwc/']; %
irun=244;  direc(irun).dir=[topdir(irun).dir 'psacw_0.35/']; %Rnc=240
irun=245;  direc(irun).dir=[topdir(irun).dir 'psacw_0.85/']; %



%array containing the machine code for the file list above - 1=PC, 2=Bezier/HPCx, 3+Newton
npess=ones(size(machineSTORE));
npess([14 16:18 23 24 27 28 31 37 39 41 43 44 46 47 48 49 50 52 53 54 55 56 57 58 62 63 ...
    64 65 67 68 72 74])=4; %4 processor 2-D runs
npess([29 30 32 33 36 38 42 48 71])=2;
npess([61])=8;

%npess([51])=30; %for 3-d runs don't need to divide TimeAv diags buy npes (see subroutine avdg)

ground_heights_store=ones(size(machineSTORE))*620;
ground_heights_store([25 34 35 36])=1000; %Miles City cases with 1000m abmsl ground height

icediag_type=ones(size(direc));

icediag_type([1 2])=2; %2 is the old type with icediag whereas 1 is the new type with icediag_nums


endt={'_1-65' '_1-88' '_1-65' '_1-88' '_1-58' '_1-88' '_1-55' '_1-60' '' '' '' '_1-88' '' '' '' '' '' '' ''};
endtALL={'' '' '_1-65' '_1-88' '_1-58' '_1-88' '_1-55' '_1-60' '' '' '' '_1-88' '' '' '' '' '' '' '' '' '' '' '' '' ''};
endtWmax={'_1-65' '_1-65' '_1-65' '_1-88' '_1-88' '' '_1-55' '_1-60' '' '' '' '_1-88' '' '' '' '' '' '' ''};

jshifts=[37:38 41:47]; %indices for 2000 km runs that need to be shifted relative to other runs
                 %because bubble was started early

if gridflag==1 & justname==0
    clear GridDan
end
 
for jdir=1:length(loadselect)   %length(direc)
    
    if justname==1
        loadflag(:)=0;
        gridflag=0;
    end
    
    j=loadselect(jdir);
    if j>12
        endt{j}='';
        endtALL{j}='';
        endtWmax{j}='';
    end
    
    direcDan(jdir).dir=direc(j).dir;
    npess2(jdir)=npess(j);
    
if justname==0        
    ground_heights(jdir)=ground_heights_store(j);    
    diagtype=icediag_type(j);
end
machine(jdir)=machineSTORE(j);

    runName(jdir).nam=run_name(j).nam;
    add_ground_heights(jdir).h=add_ground_heights_store(j).h;
    
    if strmatch(comp,'docsLap')~=1
        exdir=[direc(j).dir 'results/diags/'];
    else
        exdir=[direc(j).dir];
    end
    
    exdirSTORE(jdir).dat=exdir;
        
    
    if gridflag==1
       
		exname=[exdir 'gridDan.mat'];
		load(exname);
		GridDan(jdir)=gridDan;
    end
	
	if loadflag(1)==1
	exname=[exdir 'icediag4_5thSep2005'];
	load(exname);
	icediag4(jdir)=icediags4;
    clear icediags4;
	end
    
    if loadflag(2)==1 %microphysical process rates
	exname=[exdir 'icediag_5thSep2005'];
	load(exname);
	icediag(jdir)=icediags;
	end
	
	if loadflag(3)==1
        exname=[exdir 'maxW' endtWmax{j} '.mat'];
        load(exname);
        MaxW(jdir)=maxW;
   %     MinW(jdir)=minW;
   
       if ishift==1 & ismember(j,jshifts)
            tend=size(MaxW(jdir).w,2);
            MaxW(jdir).w(:,3:tend+2)=MaxW(jdir).w;
            MaxW(jdir).w(:,1:2)=MaxW(jdir).w(:,3:4);
        end
        
	end
    
    if loadflag(4)==1
		exname=[exdir 'vap_prctiles.mat'];
		load(exname);
		vap_prc(jdir)=vap_prctiles;
	end
    
     if loadflag(5)==1
		%exname=[exdir 'DQdehydration_10thSep2005_inc44-86'];
        %exname=[exdir 'dq_dehyd_1-' endtALL{j}];
        exname=[exdir 'dq_dehyd' endtALL{j} '.mat'];
        load(exname,'dq_vap','dq_tots');
        dq_tot(jdir)=dq_tots;
        if exist('dq_vap'); dq_vaps(jdir)=dq_vap; end
        
        if ishift==1 & ismember(j,jshifts)
            tend=size(dq_tot(jdir).d,2);
            dq_tot(jdir).d(:,3:tend+2,:)=dq_tot(jdir).d;
            dq_tot(jdir).d(:,1:2,:)=dq_tot(jdir).d(:,3:4,:);
            
            tend=size(dq_vaps(jdir).d,2);
            dq_vaps(jdir).d(:,3:tend+2,:)=dq_vaps(jdir).d;
            dq_vaps(jdir).d(:,1:2,:)=dq_vaps(jdir).d(:,3:4,:);
        end
        
        
	end
    
    if loadflag(6)==1
		exname=[exdir 'vap+tot_prcs.mat'];
		load(exname,'tot_prc','vap_prc');
        vap_prctiles(jdir)=vap_prc;
        tot_prctiles(jdir)=tot_prc;
	end
    
    if loadflag(7)==1
		exname=[exdir 'nn.mat'];
		load(exname,'n');
        nn(jdir)=n;
        load(exname,'n2');
        nn2(jdir)=n2;
        
        if ishift==1 & ismember(j,jshifts)
            tend=size(nn(jdir).n,2);
            nn(jdir).n(:,3:tend+2,:)=nn(jdir).n;
            nn(jdir).n(:,1:2,:)=nn(jdir).n(:,3:4,:);
            nn2(jdir).n(:,3:tend+2,:)=nn2(jdir).n;
            nn2(jdir).n(:,1:2,:)=nn2(jdir).n(:,3:4,:);
        end
        
        
	end
    
    if loadflag(8)==1
		exname=[exdir 'icediagsALL' endtALL{j} '.mat'];
		load(exname,'icediagALL');
        icediagsALL(jdir)=icediagALL;
        
        if ishift==1 & ismember(j,jshifts)
            tend=size(icediagsALL(jdir).i,2);
            icediagsALL(jdir).i(:,3:tend+2,:)=icediagsALL(jdir).i;
            icediagsALL(jdir).i(:,1:2,:)=icediagsALL(jdir).i(:,3:4,:);
        end
	end
    
    if loadflag(9)==1
		exname=[exdir 'icediag' endt{j}];
		load(exname,'icediags');
        icediag(jdir)=icediags;
	end
    
    if loadflag(10)==1
		exname=[exdir 'fluxes_1-65'];
        load(exname,'flux');
		fluxes(jdir)=flux;
	end
    
    if loadflag(11)==1
        if diagtype==2
            exname=[exdir 'icediag' endt{j}];
			load(exname,'icediags');
            icediag(jdir)=icediags;
        else
			exname=[exdir 'icediag_nums.mat' endt{j}];
			load(exname,'icediags_nums');
            icediag(jdir)=icediags_nums;
        end
        
         if ishift==1 & ismember(j,jshifts)
            tend=size(icediag(jdir).i,2);
            icediag(jdir).i(:,3:tend+2,:)=icediag(jdir).i;
            icediag(jdir).i(:,1:2,:)=icediag(jdir).i(:,3:4,:);
        end
        
        
	end
    
%     if loadflag(12)==1
%         exname=[exdir 'W_prctiles'];
%         load(exname,'w');
%         w_prctiles(jdir)=w;	
% 	end
    
    if loadflag(12)==1
        %exname=[exdir 'maxW_1-' endt{j}];
        clear wd wb w
        exname=[exdir 'W_prc+dist.mat'];
        load(exname,'w','wb','wd');
        
        w_prctiles(jdir)=w;
        wdist(jdir)=wd;
        wbb(jdir)=wb;
        
    end
    
    if loadflag(16)==1
        
        exname=[exdir 'dq_potemp_medacc0.05.mat'];
        
        load(exname,'dqnonSAVE','nnnonSAVE','dqpoSAVE','nnpoSAVE');
        
        dqnon(jdir)=dqnonSAVE;
        nnnon(jdir)=nnnonSAVE;
        
        dqpo(jdir)=dqpoSAVE;
        nnpo(jdir)=nnpoSAVE;
        
    end
    
     if loadflag(17)==1
         
        exname=[exdir 'simaxTimH2.mat'];
        
%        load(exname,'simaxSAVE','simeanSAVE','siminSAVE');
        load(exname,'simaxSAVE');
        
        simaxTimH(jdir)=simaxSAVE;
      %  siminTimH(jdir)=siminSAVE;
%        simean(jdir)=simeanSAVE;
        
    end
    
    if loadflag(18)==1
        
           exname=[exdir 'rhoTpert.mat'];
        
        %load(exname,'rhopertSAVE','rhopertmaxSAVE','rhopertminSAVE','rhopertSAVE2','rhopertmaxSAVE2','rhopertminSAVE2')
        
        load(exname,'rhopertSAVE','tpertSAVE');
        tpertTimH(jdir)=tpertSAVE;
        
        
    end
    
     if loadflag(19)==1
              
        exname=[exdir 'rhopertALL2.mat'];
        
        load(exname,'rhopertmaxSAVE','rhopertminSAVE','rhopertSAVE','rhopertmaxSAVE2','rhopertminSAVE2','rhopertSAVE2');
        
        rhopertmaxTimH(jdir)=rhopertmaxSAVE;
        rhopertminTimH(jdir)=rhopertminSAVE;
        rhopertTimH(jdir)=rhopertSAVE;
        
        rhopertTimHmax2(jdir)=rhopertmaxSAVE2;
        rhopertTimHmin2(jdir)=rhopertminSAVE2;
        rhopertTimH2(jdir)=rhopertSAVE2;

    end
    
    if loadflag(20)==1
        
        exname=[exdir 'LNB.mat'];
        load(exname,'me_lnb_abvSAVE','me_lnb_belSAVE','minlnbSAVE','maxlnbSAVE','lnbbinsSAVE')
        
        meanlnb_abv(jdir)=me_lnb_abvSAVE;
        meanlnb_bel(jdir)=me_lnb_belSAVE;
        
        minlnb(jdir)=minlnbSAVE;
        maxlnb(jdir)=maxlnbSAVE;
        
        lnbbins(jdir)=lnbbinsSAVE;

    end
    
    if loadflag(21)==1
        
        exname=[exdir 'LNB_bel5ppmv.mat'];
        
        load(exname,'bins_vapSAVE','bins_totSAVE','me_lnb_abv_vapSAVE','me_lnb_bel_vapSAVE','minlnb_vapSAVE','maxlnb_vapSAVE','lnbbins_pos_vapSAVE','lnbbins_neg_vapSAVE',...
            'me_lnb_abv_totSAVE','me_lnb_bel_totSAVE','minlnb_totSAVE','maxlnb_totSAVE','lnbbins_pos_totSAVE','lnbbins_neg_totSAVE');
        
        meanlnb_abv_vap(jdir)=me_lnb_abv_vapSAVE;
        meanlnb_bel_vap(jdir)=me_lnb_bel_vapSAVE;
        
        minlnb_vap(jdir)=minlnb_vapSAVE;
        maxlnb_vap(jdir)=maxlnb_vapSAVE;
        
        lnbbins_pos_vap(jdir)=lnbbins_pos_vapSAVE;
        lnbbins_neg_vap(jdir)=lnbbins_neg_vapSAVE;
        
        meanlnb_abv_tot(jdir)=me_lnb_abv_totSAVE;
        meanlnb_bel_tot(jdir)=me_lnb_bel_totSAVE;
        
        minlnb_tot(jdir)=minlnb_totSAVE;
        maxlnb_tot(jdir)=maxlnb_totSAVE;
        
        lnbbins_pos_tot(jdir)=lnbbins_pos_totSAVE;
        lnbbins_neg_tot(jdir)=lnbbins_neg_totSAVE;
        
        
        bins_vap(jdir).b=bins_vapSAVE;
        bins_tot(jdir).b=bins_totSAVE;
    end
    
	if loadflag(22)==1
        exname=[exdir 'diagsRAD.mat'];
		load(exname,'icediagsRADSAVE');
		icediagsRAD(jdir)=icediagsRADSAVE;
	end

    
    if loadflag(23)==1
		exname=[exdir 'fall_dist.mat'];
		%exname=['c:/documents and settings/g/my documents/lem/fall_dist.mat'];
        load(exname,'distSAVE');
		dist(jdir)=distSAVE;
    end
    
     if loadflag(24)==1
     	exname=[exdir 'diagsTEMP.mat'];
		load(exname,'icediagstempSAVE');
		icediagsTEMP=icediagstempSAVE(jdir);
    end
    
    if loadflag(25)==1
        exname=[exdir 'rhopert5ppmv.mat'];
        
        
        
        varlist={'rhopertTimH','rhopertTimHmax','rhopertTimHmin',...
                'rhopertTimH2','rhopertTimHmax2','rhopertTimHmin2',...
                'rho5ppmv_vap','rho5ppmv_tot',...
                'rho5ppmv_vapneg','rho5ppmv_vappos',...
                'rho5ppmv_totneg','rho5ppmv_totpos',...
                'low_prctiles','mid_prctiles','upp_prctiles','ztr_prctiles',...
                'tra'};
        
        
        savecom='load(exname';
        for i=1:length(varlist)
            savestr{i}=[varlist{i} 'SAVE'];
			savecom=[savecom ',''' savestr{i} ''''];
        end
        savecom=[savecom ');'];        
        eval(savecom);
        
        for i=1:length(varlist)
            savestr{i}=[varlist{i} 'SAVE'];
            comm=[ varlist{i} '(jdir)=' varlist{i} 'SAVE'];
            eval(comm);     
        end
        
    end

	if loadflag(26)==1
		exname=[exdir 'vaptotDists.mat'];
        varlist={'totdist','vapdist'};
%        varlist={'vapdist'};
        
        loadingCommands %saves variables in varlist (with SAVE on end)

    end
    
    
    
    if loadflag(27)==1
        exname=[exdir 'windW_1-20_short.mat'];
        varlist={'TwoD_alltim'};
        
        loadingCommands %saves variables in varlist (with SAVE on end)
    end
    
    if loadflag(29)==1
        exname=[exdir 'qprctiles.mat'];
        varlist={'q_prctiles'};
        
        loadingCommands %saves variables in varlist (with SAVE on end)
    end
    
    if loadflag(30)==1
        exname=[exdir 'dump62icediags'];
        varlist={'dump62icediags'};
        
        loadingCommands %saves variables in varlist (with SAVE on end)
    end
    
     if loadflag(31)==1
        exname=[exdir 'TRACERacc4'];
        varlist={'TRACERacc3','CONDacc3','RAINacc3','LWCacc3','rhoacc3'};
        loadingCommands %saves variables in varlist (with SAVE on end)
    end
    
    if loadflag(32)==1
        exname=[exdir 'INCgt1e8'];
        varlist={'rhoacc_alltim','INCacc_alltim','INCmaxacc_alltim','IWCacc_alltim'};
        loadingCommands %saves variables in varlist (with SAVE on end)
    end
    
    if loadflag(33)==1
        exname=[exdir 'updraught_width'];
        varlist={'width'};
        loadingCommands %saves variables in varlist (with SAVE on end)
    end
    
     if loadflag(34)==1
         switch comp
         case 'lacieLap'
            exname=['c:/documents and settings/g/my documents/balloon_diags'];
         case 'uni'
            exname=['c:/documents and settings/login/my documents/hibiscus/baloondata/balloon_diags'];
         end
        varlist={'dmi','sdla','saw'};
        loadingCommands %saves variables in varlist (with SAVE on end)
    end
    
     if loadflag(35)==1
        exname=[exdir 'TRACERmax'];
        varlist={'TRACERmax'};
        loadingCommands %saves variables in varlist (with SAVE on end)
    end
    
    if loadflag(36)==1
        exname=[exdir 'rho_mean.mat'];
        varlist={'rho_prof'};
        loadingCommands %saves variables in varlist (with SAVE on end)
    end
    
    if loadflag(37)==1
        exname=[exdir 'radar_10dbz.mat'];
        varlist={'n10dbz'};
        loadingCommands %saves variables in varlist (with SAVE on end)
    end
    
    if loadflag(38)==1 %ARM Darwin echo tops
        exname=[exdir 'radar_echotops.mat'];
        varlist={'ntop'};
        loadingCommands %saves variables in varlist (with SAVE on end)
    end
    
    if loadflag(39)==1 %temperature perturbations
        exname=[exdir 'Tperts_all.mat'];
        varlist={'tpertTimH_full'};
        loadingCommands %saves variables in varlist (with SAVE on end)
    end
    
    if loadflag(40)==1
         exname=[exdir 'turbq_flux.mat'];
        varlist={'vappos_up',
                'vapneg_up',
                'vappos_d',
                'vapneg_d',
                'totpos_up',
                'totneg_up',
                'totpos_d',  %totpod_d and vappos_d might be missing
                'totneg_d' };
                    
        loadingCommands %saves variables in varlist (with SAVE on end)
    end
    
    if loadflag(41)==1 %temperature perturbations
        exname=[exdir 'vapfull.mat'];
        varlist={'vapfull'};
        loadingCommands %saves variables in varlist (with SAVE on end)
    end
    
    if loadflag(42)==1
        exname=[exdir 'fields_for_supersat_hslices2.mat'];
        varlist={'pressure_hslice','vap_hslice','potemp_hslice'};
        loadingCommands %saves variables in varlist (with SAVE on end)
    end 
    
	if loadflag(43)==1
        exname=[exdir 'ice_hslices_justice.mat'];
        varlist={'iceMR_hslice','iceNC_hslice'}; %,'snowMR_hslice','snowNC_hslice'};
        loadingCommands %saves variables in varlist (with SAVE on end)
    end 
    
     if loadflag(45)==1
        exname=[exdir 'totice_44.mat'];
        varlist={'totice_44'}; %,'snowMR_hslice','snowNC_hslice'};
        loadingCommands %saves variables in varlist (with SAVE on end)
    end 
    
    if loadflag(46)==1
        exname=[exdir 'dT.mat'];
        varlist={'dT_conv','dT_nonconv','dT_bubble'}; 
        loadingCommands %saves variables in varlist (with SAVE on end)
    end 
    if loadflag(47)==1
        exname=[exdir 'tperts2.mat'];
        varlist={'tpertTimH_full','tpert_max','vappert_max','vappertTimH_full','capes'}; 
        loadingCommands %saves variables in varlist (with SAVE on end)
    end     
    
 
end    %for jdir=1:length(loadselect) 
'loadVapData finished'

%