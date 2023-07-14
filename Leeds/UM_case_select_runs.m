function [vars_out]=UM_case_select_runs(UM_cases,iadd_umid_label)

if ~exist('iadd_umid_label')
    iadd_umid_label=1;
end

clear dirUM  %Clearing this, so that it can be set as either a cell (for multiple specifications of 
% the root directory for the UM files) or a string.

for idat=1:99
    flag{idat}='';  %The default is netCDF
    i_aerosol_processing_multi{idat}=0;
    append_str_timser{idat}='';
end

%UM_cases = '12th Nov case, as of May 2016';


switch UM_cases
    case '12th Nov case, as of May 2016'
        %% Newest runs May 2016


        dirUM='/home/disk/eos8/d.grosvenor/UM/26thOct_POC/';
        dirUM='/home/disk/eos8/d.grosvenor/UM/12Nov2008_Boutle/';

        idat=1;
        fileUM{idat} = '/xmmz-u/xmmzu_VAR_NAME_.pp.nc'; labs_UM(idat).l = 'CASIM-Ndvar';  flag{idat} = 'load_mat'; fileUM_rho{idat} = '/xmmz-u/xmmzu_rho_.pp.nc';pole_lat=70; pole_lon=284;
        line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.4 0.4 0.4]; marker_styleUM(idat).m='d'; idat=idat+1;
        fileUM{idat} = '/xmmz-v/xmmzv_VAR_NAME_.pp.nc'; labs_UM(idat).l ='CASIM-Ndvar-0.1';  flag{idat} = 'load_mat'; fileUM_rho{idat} = '/xmmz-v/xmmzv_rho_.pp.nc'; pole_lat=70; pole_lon=284;
        line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0 0 0]; marker_styleUM(idat).m='v'; idat=idat+1;
        fileUM{idat} = '/xmmz-w/xmmzw_VAR_NAME_.pp.nc'; labs_UM(idat).l ='CASIM-Ndvar-10';  flag{idat} = 'load_mat'; fileUM_rho{idat} = '/xmmz-v/xmmzv_rho_.pp.nc'; pole_lat=70; pole_lon=284;
        line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.8 0.4 0]; marker_styleUM(idat).m='^'; idat=idat+1;
        %fileUM{idat} = '/xlhg-v/xlhgv_VAR_NAME_.pp.nc'; labs_UM(idat).l = 'CASIM-Ndvar-0.025';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
        %        line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0 0.8 0]; marker_styleUM(idat).m='v'; idat=idat+1;
        fileUM{idat} = '/xmmz-x/xmmzx_VAR_NAME_.pp.nc'; labs_UM(idat).l = 'CASIM-Ndvar-0.025';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
        line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0 0.8 0]; marker_styleUM(idat).m='o'; idat=idat+1;
        fileUM{idat} = '/xmmz-n/xmmzn_VAR_NAME_.pp.nc'; labs_UM(idat).l = 'Old-mphys';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
        line_patternUM(idat).p= '--';  line_colourUM(idat).c=[1 0 0]; marker_styleUM(idat).m='s'; idat=idat+1;
        %fileUM{idat} = '/xmmz-l/xmmzl_VAR_NAME_.pp.nc'; labs_UM(idat).l ='CASIM-Ndvar-RHcrit0.999';  flag{idat} = 'load_mat'; fileUM_rho{idat} = '/xmmz-v/xmmzv_rho_.pp.nc'; pole_lat=70; pole_lon=284;
        %    line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0 0 0]; marker_styleUM(idat).m='v'; idat=idat+1;
    
    case '12th Nov case, as of May 2016 adhoc'
        %% Newest runs May 2016


        dirUM='/home/disk/eos8/d.grosvenor/UM/26thOct_POC/';
        dirUM='/home/disk/eos8/d.grosvenor/UM/12Nov2008_Boutle/';

        idat=1;
        fileUM{idat} = '/xmmz-u/xmmzu_VAR_NAME_.pp.nc'; labs_UM(idat).l = 'CASIM-Ndvar';  flag{idat} = 'load_mat'; fileUM_rho{idat} = '/xmmz-u/xmmzu_rho_.pp.nc';pole_lat=70; pole_lon=284;
        line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.4 0.4 0.4]; marker_styleUM(idat).m='d'; idat=idat+1;
%        fileUM{idat} = '/xmmz-v/xmmzv_VAR_NAME_.pp.nc'; labs_UM(idat).l ='CASIM-Ndvar-0.1';  flag{idat} = 'load_mat'; fileUM_rho{idat} = '/xmmz-v/xmmzv_rho_.pp.nc'; pole_lat=70; pole_lon=284;
%        line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0 0 0]; marker_styleUM(idat).m='v'; idat=idat+1;
%        fileUM{idat} = '/xmmz-w/xmmzw_VAR_NAME_.pp.nc'; labs_UM(idat).l ='CASIM-Ndvar-10';  flag{idat} = 'load_mat'; fileUM_rho{idat} = '/xmmz-v/xmmzv_rho_.pp.nc'; pole_lat=70; pole_lon=284;
%        line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.8 0.4 0]; marker_styleUM(idat).m='^'; idat=idat+1;
%         fileUM{idat} = '/xlhg-v/xlhgv_VAR_NAME_.pp.nc'; labs_UM(idat).l = 'CASIM-Ndvar-0.025';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
%         line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0 0.8 0]; marker_styleUM(idat).m='v'; idat=idat+1;
%        fileUM{idat} = '/xmmz-x/xmmzx_VAR_NAME_.pp.nc'; labs_UM(idat).l = 'CASIM-Ndvar-0.025';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
%        line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0 0.8 0]; marker_styleUM(idat).m='o'; idat=idat+1;
%        fileUM{idat} = '/xmmz-n/xmmzn_VAR_NAME_.pp.nc'; labs_UM(idat).l = 'Old-mphys';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
%        line_patternUM(idat).p= '--';  line_colourUM(idat).c=[1 0 0]; marker_styleUM(idat).m='s'; idat=idat+1;
        %fileUM{idat} = '/xmmz-l/xmmzl_VAR_NAME_.pp.nc'; labs_UM(idat).l ='CASIM-Ndvar-RHcrit0.999';  flag{idat} = 'load_mat'; fileUM_rho{idat} = '/xmmz-v/xmmzv_rho_.pp.nc'; pole_lat=70; pole_lon=284;
        %    line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0 0 0]; marker_styleUM(idat).m='v'; idat=idat+1;
%        fileUM{idat} = '/xmmz-p/xmmzp_VAR_NAME_.pp.nc'; labs_UM(idat).l = 'CASIM-Ndvar-sed';  flag{idat} = 'load_mat'; fileUM_rho{idat} = '';pole_lat=70; pole_lon=284;
%        line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.4 0.4 0.4]; marker_styleUM(idat).m='d'; idat=idat+1;    
%        fileUM{idat} = '/xmmz-q/xmmzq_VAR_NAME_.pp.nc'; labs_UM(idat).l = 'CASIM-Ndvar-10-sed';  flag{idat} = 'load_mat'; fileUM_rho{idat} = '';pole_lat=70; pole_lon=284;
%        line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.4 0.4 0.4]; marker_styleUM(idat).m='d'; idat=idat+1;                   
%        fileUM{idat} = '/xmmz-j/xmmzj_VAR_NAME_.pp.nc'; labs_UM(idat).l = 'CASIM-fixed500_CloudScheme_OFF';  flag{idat} = 'load_mat'; fileUM_rho{idat} = '';pole_lat=70; pole_lon=284;
%        line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.4 0.4 0.4]; marker_styleUM(idat).m='d'; idat=idat+1;                           

    case '12th Nov case, as of May 2016 adhoc eos10'
        
        dirUM='/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/';        

        idat=1;
        
%        fileUM{idat} = '/xmmz-t/xmmzt_VAR_NAME_.pp.nc'; labs_UM(idat).l ='CASIM-Ndvar-lmrphys';  flag{idat} = 'load_mat'; fileUM_rho{idat} = '/xmmz-t/xmmzt_rho_.pp.nc'; pole_lat=70; pole_lon=284;
%        line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.8 0.4 0]; marker_styleUM(idat).m='^'; idat=idat+1;

        fileUM{idat} = 'u-ag688/u-ag688_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc'; labs_UM(idat).l ='CASIMv10.3-Ndvar';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
        line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.8 0.4 0]; marker_styleUM(idat).m='^'; idat=idat+1;

    case '12th Nov case, as of May 2016 adhoc multi-dirUM'
        
        dirUM2{1}='/home/disk/eos8/d.grosvenor/UM/12Nov2008_Boutle/';                
        dirUM2{2}='/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/';        

        idat=1;
        
%        fileUM{idat} = '/xmmz-t/xmmzt_VAR_NAME_.pp.nc'; labs_UM(idat).l ='CASIM-Ndvar-lmrphys';  flag{idat} = 'load_mat'; fileUM_rho{idat} = '/xmmz-t/xmmzt_rho_.pp.nc'; pole_lat=70; pole_lon=284;
%        line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.8 0.4 0]; marker_styleUM(idat).m='^'; idat=idat+1;

% Rain on runs
%         fileUM{idat} = 'u-ah963/u-ah963_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc'; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='CASIMv10.3-Ndvar';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
%         line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.8 0.4 0]; marker_styleUM(idat).m='^'; idat=idat+1;
%         
%         fileUM{idat} = 'u-ag688/u-ag688_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc'; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='CASIMv10.3-Ndvar-proc';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
%         line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.8 0.4 0]; marker_styleUM(idat).m='^'; idat=idat+1;  
%         
%         fileUM{idat} = 'u-ah971/u-ah971_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc'; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='CASIMv10.3-Ndvar-cloudSchemeOFF';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
%         line_patternUM(idat).p= '--';  line_colourUM(idat).c=[1 0.0 0.0]; marker_styleUM(idat).m='^'; idat=idat+1;
%         
%         um_str{idat}='u-ah964'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='CASIMv10.3-Ndvar-proc-cloudSchemeOFF';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
%         line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.4 0.4 0.4]; marker_styleUM(idat).m='^'; idat=idat+1;        

%Rain off runs        
%       um_str{idat}='u-ai120'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='CASIMv10.3-Ndvar-no-auto';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
%       line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.8 0.4 0]; marker_styleUM(idat).m='^'; idat=idat+1;
        
%       um_str{idat}='u-ai117'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='CASIMv10.3-Ndvar-no-auto-proc';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
%       line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.8 0.4 0]; marker_styleUM(idat).m='^'; idat=idat+1;
  
%        um_str{idat}='u-ai121'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='CASIMv10.3-Ndvar-no-auto-proc-cloudSchemeOFF';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
%        line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.4 0.4 0.4]; marker_styleUM(idat).m='^'; idat=idat+1;
%         
%          um_str{idat}='u-ai123'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='CASIMv10.3-Ndvar-no-auto-cloudSchemeOFF';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
%          line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.4 0.4 0.4]; marker_styleUM(idat).m='^'; idat=idat+1;             
% % 
%        um_str{idat}='u-ai397'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='CASIMv10.3-Ndvar-proc-cloudSchemeOFF-no-auto-surf-source';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
%        line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.4 0.4 0.4]; marker_styleUM(idat).m='^'; idat=idat+1;        

%          um_str{idat}='u-ai417'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='CASIMv10.3-Ndvar-Annette';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
%          line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.4 0.8 0.4]; marker_styleUM(idat).m='^'; idat=idat+1;             
%          
%          um_str{idat}='u-ai527'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='CASIMv10.3-Ndvar-Annette-BL-mixing-OFF';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
%          line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.8 0.8 0.8]; marker_styleUM(idat).m='^'; idat=idat+1;
% 
%          um_str{idat}='u-ai528'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='CASIMv10.3-Ndvar-Annette-Mass-cons';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
%          line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.8 0.4 0.4]; marker_styleUM(idat).m='s'; idat=idat+1;
% 
%          um_str{idat}='u-ai529'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='CASIMv10.3-Ndvar-Process-level1-sed-changes';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
%          line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.4 0.4 0.8]; marker_styleUM(idat).m='s'; idat=idat+1;
% 
%          um_str{idat}='u-ai531'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='CASIMv10.3-Ndvar-l-mr-physics1-true';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
%          line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.8 0.8 0.4]; marker_styleUM(idat).m='^'; idat=idat+1;

         %um_str{idat}='u-ai535'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='CASIMv10.3-Ndvar-sed-changes';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
         %line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.4 0.4 0.8]; marker_styleUM(idat).m='o'; idat=idat+1;

         

         %um_str{idat}='u-aj091'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='CASIMv10.3-Ndvar-Processing_No_mix_tr_bl';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
         %line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.8 0.8 0.4]; marker_styleUM(idat).m='^'; idat=idat+1;

         %um_str{idat}='u-aj097'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='CASIMv10.3-Ndvar-No_mix_tr_bl';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
         %line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.4 0.4 0.8]; marker_styleUM(idat).m='o'; idat=idat+1;
%           
        
%          um_str{idat}='u-aj470'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='CASIMv10.3-Ndvar-Processing_No_mix_tr_bl';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
%          line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.8 0.8 0.4]; marker_styleUM(idat).m='^'; idat=idat+1;
% 
%          um_str{idat}='u-aj471'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='CASIMv10.3-Ndvar-No_mix_tr_bl';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
%          line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.4 0.4 0.8]; marker_styleUM(idat).m='o'; idat=idat+1;
% 
%          um_str{idat}='u-aj472'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='CASIMv10.3-Ndvar-Processing_No_mix_tr_bl';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
%          line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.8 0.8 0.4]; marker_styleUM(idat).m='^'; idat=idat+1;
% 
%          um_str{idat}='u-aj475'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='CASIMv10.3-Ndvar-No_mix_tr_bl';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
%          line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.4 0.4 0.8]; marker_styleUM(idat).m='o'; idat=idat+1;
%          
%          um_str{idat}='u-aj476'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='CASIMv10.3-Ndvar-Processing_No_mix_tr_bl';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
%          line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.8 0.8 0.4]; marker_styleUM(idat).m='^'; idat=idat+1;


        um_str{idat}='xmmz-u'; fileUM{idat} = [um_str{idat} '/xmmzu_VAR_NAME_.pp.nc']; dirUM{idat}=dirUM2{1}; labs_UM(idat).l = 'CASIMv8.5-Ndvar-210cm3';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
        line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.4 0.4 0.4]; marker_styleUM(idat).m='d'; idat=idat+1;
     



case '12th Nov case, as of May 2016 processing runs multi-dirUM'
        
        dirUM2{1}='/home/disk/eos8/d.grosvenor/UM/12Nov2008_Boutle/';                
        dirUM2{2}='/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/';        

        idat=1;
        
%        fileUM{idat} = '/xmmz-t/xmmzt_VAR_NAME_.pp.nc'; labs_UM(idat).l ='CASIM-Ndvar-lmrphys';  flag{idat} = 'load_mat'; fileUM_rho{idat} = '/xmmz-t/xmmzt_rho_.pp.nc'; pole_lat=70; pole_lon=284;
%        line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.8 0.4 0]; marker_styleUM(idat).m='^'; idat=idat+1;

        fileUM{idat} = 'u-ah963/u-ah963_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc'; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='CASIMv10.3-Ndvar';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
        line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.8 0.4 0]; marker_styleUM(idat).m='^'; idat=idat+1;
        
        fileUM{idat} = 'u-ag688/u-ag688_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc'; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='CASIMv10.3-Ndvar-proc';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
        line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.8 0.4 0]; marker_styleUM(idat).m='^'; idat=idat+1;
  
        fileUM{idat} = 'u-ah964/u-ah964_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc'; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='CASIMv10.3-Ndvar-proc-cloudSchemeOFF';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
        line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.4 0.4 0.4]; marker_styleUM(idat).m='^'; idat=idat+1;
        
        fileUM{idat} = 'u-ah971/u-ah971_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc'; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='CASIMv10.3-Ndvar-cloudSchemeOFF';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
        line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.4 0.4 0.4]; marker_styleUM(idat).m='^'; idat=idat+1;
               
        
%        um_str{idat}='u-ai397'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='CASIMv10.3-Ndvar-proc-no-rain-surf-source';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
%        line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.4 0.4 0.4]; marker_styleUM(idat).m='s'; idat=idat+1;                
        
case '12th Nov case, as of May 2016 processing runs PLOTS multi-dirUM'
        
        dirUM2{1}='/home/disk/eos8/d.grosvenor/UM/12Nov2008_Boutle/';                
        dirUM2{2}='/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/';        

        idat=1;


        um_str{idat}='xmmz-u'; fileUM{idat} = [um_str{idat} '/xmmzu_VAR_NAME_.pp.nc']; dirUM{idat}=dirUM2{1}; labs_UM(idat).l = 'CASIMv8.5-Ndvar-210cm3';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
        line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.4 0.4 0.4]; marker_styleUM(idat).m='d'; idat=idat+1;
        
        um_str{idat}='u-ah963'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='CASIMv10.3-Ndvar-600cm3';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
        line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.8 0.4 0]; marker_styleUM(idat).m='^'; idat=idat+1;
        
        um_str{idat}='u-aj091'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='CASIMv10.3-Ndvar-210cm3-No-mix-tr-bl';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
        line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.8 0.8 0.4]; marker_styleUM(idat).m='^'; idat=idat+1;
        
         um_str{idat}='u-aj097'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='CASIMv10.3-Ndvar-210cm3-Processing-No-mix-tr-bl';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
         line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.4 0.4 0.8]; marker_styleUM(idat).m='o'; idat=idat+1;

         um_str{idat}='u-ai527'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='CASIMv10.3-Ndvar-600cm3-Annette-BL-mixing-OFF';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
         line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.8 0.8 0.8]; marker_styleUM(idat).m='^'; idat=idat+1;
        
         fileUM{idat} = 'u-ag688/u-ag688_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc'; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='CASIMv10.3-Ndvar-600cm3-proc';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
         line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.8 0.4 0]; marker_styleUM(idat).m='^'; idat=idat+1;
%         

    case '12th Nov case, as of Feb 2017 processing runs PLOTS multi-dirUM'
        
        %dirUM2{1}='/home/disk/eos8/d.grosvenor/UM/12Nov2008_Boutle/'; %Not used here
        dirUM2{2}='/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/';        

        idat=1;
% %         um_str{idat}='u-aj091'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='CASIMv10.3-Ndvar-No-mix-tr-bl-210cm3';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
% %         line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.8 0.8 0.4]; marker_styleUM(idat).m='^'; idat=idat+1;
% %         
% %          um_str{idat}='u-aj097'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='CASIMv10.3-Ndvar-Processing-No-mix-tr-bl-210cm3';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
% %          line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.4 0.4 0.8]; marker_styleUM(idat).m='o'; idat=idat+1;
% %         
% %          um_str{idat}='u-aj470'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='CASIMv10.3-Ndvar-Processing-No-mix-tr-bl-600cm3';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
% %          line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.8 0.4 0.4]; marker_styleUM(idat).m='^'; idat=idat+1;
% % 
%            um_str{idat}='u-aj471'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='CASIMv10.3-Ndvar-Processing-No-mix-tr-bl-210cm3-CS-OFF';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
%            line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.4 0.4 0.4]; marker_styleUM(idat).m='o'; idat=idat+1;
% % % 
%            um_str{idat}='u-aj472'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='CASIMv10.3-Ndvar-No-mix-tr-bl-210cm3-CS-OFF';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
%            line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.8 0.8 0.8]; marker_styleUM(idat).m='^'; idat=idat+1;
% % % 
% %          um_str{idat}='u-aj475'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='CASIMv10.3-Ndvar-Processing-No-mix-tr-bl-210cm3-rain-OFF';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
% %          line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.4 0.8 0.8]; marker_styleUM(idat).m='o'; idat=idat+1;
% %          
% %          um_str{idat}='u-aj476'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='CASIMv10.3-Ndvar-No-mix-tr-bl-210cm3-rain-OFF';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
% %          line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.8 0.4 0.8]; marker_styleUM(idat).m='^'; idat=idat+1;
% 
%           um_str{idat}='u-al076'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='CASIMv10.3-Ndvar-Processing-No-mix-tr-bl-210cm3-CS-OFF-rain-OFF-Drop-sed-OFF-No-coarse';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
%           line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.4 0.8 0.8]; marker_styleUM(idat).m='o'; idat=idat+1;
% %          
%          um_str{idat}='u-al119'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='CASIMv10.3-Ndvar-Processing-No-mix-tr-bl-210cm3-CS-OFF-rain-OFF-BUG-FIX-N-Drop-sed-OFF-No-coarse';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
%          line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.8 0.4 0.8]; marker_styleUM(idat).m='^'; idat=idat+1;

% um_str{idat}='u-al527'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='Controlx0.025-RainOFF-CloudSchemeON';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
% line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.0 0.8 0.0]; marker_styleUM(idat).m='^'; idat=idat+1;
% 
%um_str{idat}='u-al136'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='Control-Processing-tests';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284; %'Control-Processing-No-mix-tr-bl-CS-OFF-rain-OFF-BUG-FIX-MASS-N-Drop-sed-OFF-No-coarse'
%line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.8 0.8 0.8]; marker_styleUM(idat).m='^'; idat=idat+1;
% 
%

um_str{idat}='u-am784'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='Constant flux 1e6';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.4 0.4 0.4]; marker_styleUM(idat).m='^'; idat=idat+1;

%um_str{idat}='u-al522'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='Control-Processing';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
%line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.4 0.4 0.4]; marker_styleUM(idat).m='^'; idat=idat+1;

% um_str{idat}='u-al523'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='Control';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
% line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.4 0.4 0.4]; marker_styleUM(idat).m='^'; idat=idat+1;
% 
% um_str{idat}='u-al524'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='Controlx0.1-Processing';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
% line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.0 0.0 0.0]; marker_styleUM(idat).m='^'; idat=idat+1;
% 
% um_str{idat}='u-al525'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='Controlx0.1';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
% line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.0 0.0 0.0]; marker_styleUM(idat).m='^'; idat=idat+1;
% 
% um_str{idat}='u-al526'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='Controlx0.025-Processing';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
% line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.0 0.8 0.0]; marker_styleUM(idat).m='^'; idat=idat+1;
% 
% um_str{idat}='u-al661'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='Controlx0.025';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
% line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.0 0.8 0.0]; marker_styleUM(idat).m='^'; idat=idat+1;
% 
% 
% um_str{idat}='u-al539'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='Controlx10-Processing';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
% line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.8 0.4 0.0]; marker_styleUM(idat).m='^'; idat=idat+1;
% 
% um_str{idat}='u-al541'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='Controlx10';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
% line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.8 0.4 0.0]; marker_styleUM(idat).m='^'; idat=idat+1;


    case '12th Nov case, as of Feb 2017 processing runs PLOTS multi-dirUM processing ON'
        
        %dirUM2{1}='/home/disk/eos8/d.grosvenor/UM/12Nov2008_Boutle/'; %Not used here
        dirUM2{2}='/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/';        

        idat=1;

% 
um_str{idat}='u-al522'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='Control-Processing';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.4 0.4 0.4]; marker_styleUM(idat).m='^'; idat=idat+1;

%um_str{idat}='u-al523'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='Control';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
%line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.4 0.4 0.4]; marker_styleUM(idat).m='^'; idat=idat+1;

um_str{idat}='u-al524'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='Controlx0.1-Processing';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.0 0.0 0.0]; marker_styleUM(idat).m='^'; idat=idat+1;

%um_str{idat}='u-al525'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='Controlx0.1';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
%line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.0 0.0 0.0]; marker_styleUM(idat).m='^'; idat=idat+1;

um_str{idat}='u-al526'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='Controlx0.025-Processing';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.0 0.8 0.0]; marker_styleUM(idat).m='^'; idat=idat+1;

%um_str{idat}='u-al661'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='Controlx0.025';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
%line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.0 0.8 0.0]; marker_styleUM(idat).m='^'; idat=idat+1;


um_str{idat}='u-al539'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='Controlx10-Processing';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.8 0.4 0.0]; marker_styleUM(idat).m='^'; idat=idat+1;

%um_str{idat}='u-al541'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='Controlx10';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
%line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.8 0.4 0.0]; marker_styleUM(idat).m='^'; idat=idat+1;


case '12th Nov case, as of Feb 2017 processing runs PLOTS multi-dirUM Processing OFF'
        
        %dirUM2{1}='/home/disk/eos8/d.grosvenor/UM/12Nov2008_Boutle/'; %Not used here
        dirUM2{2}='/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/';        

        idat=1;
%um_str{idat}='u-al136'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='Control-Processing-tests';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284; %'Control-Processing-No-mix-tr-bl-CS-OFF-rain-OFF-BUG-FIX-MASS-N-Drop-sed-OFF-No-coarse'
%line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.8 0.8 0.8]; marker_styleUM(idat).m='^'; idat=idat+1;
% 
% 
%um_str{idat}='u-al522'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='Control-Processing';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
%line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.4 0.4 0.4]; marker_styleUM(idat).m='^'; idat=idat+1;

um_str{idat}='u-al523'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='Control';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.4 0.4 0.4]; marker_styleUM(idat).m='^'; idat=idat+1;

%um_str{idat}='u-al524'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='Controlx0.1-Processing';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
%line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.0 0.0 0.0]; marker_styleUM(idat).m='^'; idat=idat+1;

um_str{idat}='u-al525'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='Controlx0.1';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.0 0.0 0.0]; marker_styleUM(idat).m='^'; idat=idat+1;

%um_str{idat}='u-al526'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='Controlx0.025-Processing';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
%line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.0 0.8 0.0]; marker_styleUM(idat).m='^'; idat=idat+1;

um_str{idat}='u-al661'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='Controlx0.025';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.0 0.8 0.0]; marker_styleUM(idat).m='^'; idat=idat+1;


%um_str{idat}='u-al539'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='Controlx10-Processing';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
%line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.8 0.4 0.0]; marker_styleUM(idat).m='^'; idat=idat+1;

um_str{idat}='u-al541'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='Controlx10';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.8 0.4 0.0]; marker_styleUM(idat).m='^'; idat=idat+1;



    case '12th Nov case, as of May 2017 processing runs surface fluxes, delaying processing, etc.'
        
        %dirUM2{1}='/home/disk/eos8/d.grosvenor/UM/12Nov2008_Boutle/'; %Not used here
        dirUM2{2}='/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/';        

        idat=1;
        %
        
        delaying_processing_case = 'per m3 namelist';        
        delaying_processing_case = 'per kg namelist';
        
        switch delaying_processing_case
            case 'per m3 namelist'
       
% %         
%   um_str{idat}='u-al523'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='Control';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
%   line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.4 0.4 0.4]; marker_styleUM(idat).m='v'; i_aerosol_processing_multi{idat}=0; idat=idat+1;
% % % 
% % % um_str{idat}='u-al522'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='Control-Processing';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
% % %  line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.4 0.4 0.4]; marker_styleUM(idat).m='^'; i_aerosol_processing_multi{idat}=1; idat=idat+1;
% % % % 
% % % % %   - N.B. all of these are based on "Control-Processing" run above
% % % % %   (u-al522) - i.e. contorl aeroosl with processing on, but with various
% % % % %   changes (e.g. surface fluxes of aerosol, etc.)
% % % % 
% % % % um_str{idat}='u-am503'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='Replenish-aerosol-every-timestep';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
% % % % line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.0 0.8 0.0]; marker_styleUM(idat).m='o'; i_aerosol_processing_multi{idat}=1; idat=idat+1;
% % % % 
% % % % um_str{idat}='u-am516'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='No-rain-first-6hrs';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
% % % % line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.0 0.8 0.0]; marker_styleUM(idat).m='s'; i_aerosol_processing_multi{idat}=1; idat=idat+1;
% % % % 
% % % %  um_str{idat}='u-am784'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='Constant flux 1e6';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
% % % %  line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.8 0.4 0.0]; marker_styleUM(idat).m='v'; i_aerosol_processing_multi{idat}=1; idat=idat+1;
% % % %  
% % % %  um_str{idat}='u-am593'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='Constant flux 1e7';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
% % % %  line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.8 0.4 0.0]; marker_styleUM(idat).m='^'; i_aerosol_processing_multi{idat}=1; idat=idat+1;
% % % % 
% % % %  um_str{idat}='u-am595'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='1e7 dropping to zero';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
% % % %  line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.0 0.0 0.0]; marker_styleUM(idat).m='o'; i_aerosol_processing_multi{idat}=1; idat=idat+1;
% % % % % 
% % % %  um_str{idat}='u-am786'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='5e6 dropping to zero';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
% % % %  line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.0 0.0 0.0]; marker_styleUM(idat).m='s'; i_aerosol_processing_multi{idat}=1; idat=idat+1;
% % % % % 
% % % %  um_str{idat}='u-am889'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='5e6 dropping to zero, Aitken+Coarse';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
% % % %  line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.8 0.0 0.0]; marker_styleUM(idat).m='v'; i_aerosol_processing_multi{idat}=1; idat=idat+1;
% % % % % 
% % % % % 
% % % % % um_str{idat}='u-al539'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='Controlx10-Processing';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
% % % % % line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.8 0.4 0.0]; marker_styleUM(idat).m='^'; i_aerosol_processing_multi{idat}=1; idat=idat+1;
% % % % % 
% % % % % um_str{idat}='u-al541'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='Controlx10';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
% % % % % line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.8 0.4 0.0]; marker_styleUM(idat).m='^'; i_aerosol_processing_multi{idat}=1; idat=idat+1;
% % % % 
%  um_str{idat}='u-an149'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='Rain OFF';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
%  line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.0 0.8 0.0]; marker_styleUM(idat).m='o'; i_aerosol_processing_multi{idat}=1; idat=idat+1;
% % 
% % % um_str{idat}='u-an171'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='\mu_{rain}=1.5';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
% % % line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.0 0.8 0.0]; marker_styleUM(idat).m='s'; i_aerosol_processing_multi{idat}=1; idat=idat+1;
% % % 
% % %  um_str{idat}='u-an172'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='\mu_{rain}=3.5';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
% % %  line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.8 0.4 0.0]; marker_styleUM(idat).m='v'; i_aerosol_processing_multi{idat}=1; idat=idat+1;
% % %  
% % % % um_str{idat}='u-an173'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='Constant flux 1e7';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
% % % % line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.8 0.4 0.0]; marker_styleUM(idat).m='^'; i_aerosol_processing_multi{idat}=1; idat=idat+1;
% % % 
% % %  um_str{idat}='u-an174'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='\sigma_{aerosol}=2.5';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
% % %  line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.0 0.0 0.0]; marker_styleUM(idat).m='o'; i_aerosol_processing_multi{idat}=1; idat=idat+1;
% % % 
% % % um_str{idat}='u-an175'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='\sigma_{aerosol}=0.5';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
% % %  line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.0 0.0 0.0]; marker_styleUM(idat).m='s'; i_aerosol_processing_multi{idat}=1; idat=idat+1;
% % % % 
% % % % um_str{idat}='u-an295'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='Constant flux 1e6, tracer mixing';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
% % % % line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.8 0.0 0.0]; marker_styleUM(idat).m='v'; i_aerosol_processing_multi{idat}=1; idat=idat+1;
% % % %  
% % % %   um_str{idat}='u-an308'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='Constant flux 1e7, tracer mixing';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
% % % %  line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.8 0.0 0.0]; marker_styleUM(idat).m='v'; i_aerosol_processing_multi{idat}=1; idat=idat+1;
% % % % % %          
% % 
% %  um_str{idat}='u-ao226'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='Midway transfer';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
% %  line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.0 0.0 0.0]; marker_styleUM(idat).m='o'; i_aerosol_processing_multi{idat}=1; idat=idat+1;
% % 
% % um_str{idat}='u-ao227'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='No coarse';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
% %  line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.0 0.0 0.0]; marker_styleUM(idat).m='s'; i_aerosol_processing_multi{idat}=1; idat=idat+1;
% % 
% % um_str{idat}='u-ao229'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='Replenishment, no coarse';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
% % line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.8 0.0 0.0]; marker_styleUM(idat).m='v'; i_aerosol_processing_multi{idat}=1; idat=idat+1;
% %  
% %   um_str{idat}='u-ao260'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='Control + bug fix';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
% %  line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.8 0.8 0.0]; marker_styleUM(idat).m='^'; i_aerosol_processing_multi{idat}=1; idat=idat+1;
%   
% %   um_str{idat}='u-ao722'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='Bug, dm1=dm';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
% %  line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.8 0.8 0.8]; marker_styleUM(idat).m='^'; i_aerosol_processing_multi{idat}=1; idat=idat+1;
%   
%  %    um_str{idat}='u-ao771'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='Bug fix, dm1=dm';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
%  % line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.0 0.0 0.8]; marker_styleUM(idat).m='^'; i_aerosol_processing_multi{idat}=1; idat=idat+1;
%   
% %     um_str{idat}='u-ao772'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='Bug fix, dm1=dm, aero thresh';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
% %  line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.8 0.0 0.8]; marker_styleUM(idat).m='o'; i_aerosol_processing_multi{idat}=1; idat=idat+1;
%   
% %      um_str{idat}='u-ao880'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='Bug fix, dm=0, aero thresh';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
% %   line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.0 0.8 0.0]; marker_styleUM(idat).m='v'; i_aerosol_processing_multi{idat}=1; idat=idat+1;
% %   
% %      um_str{idat}='u-ao883'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='Bypass which mode, aero thresh';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
% %   line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.8 0.8 0.0]; marker_styleUM(idat).m='s'; i_aerosol_processing_multi{idat}=1; idat=idat+1;
%   
% %     um_str{idat}='u-ao889'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='Bug fix, dm1=dm, aero thresh, truncate neg';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
% %  line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.8 0.0 0.0];
% %  marker_styleUM(idat).m='^'; i_aerosol_processing_multi{idat}=1; idat=idat+1;
% 
%   
% %       um_str{idat}='u-ap287'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='Bug fix, neg dm fix';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
% %   line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.8 0.0 0.8]; marker_styleUM(idat).m='o'; i_aerosol_processing_multi{idat}=1; append_str_timser{idat}='.txt'; idat=idat+1;
% %   
% %      um_str{idat}='u-ap288'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='Bypass which mode, bug fix, neg dm fix';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
% %   line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.0 0.8 0.0]; marker_styleUM(idat).m='^'; i_aerosol_processing_multi{idat}=1; append_str_timser{idat}='.txt'; idat=idat+1;
% 
%      um_str{idat}='u-an966'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='Annette, her switches, rain OFF';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
%   line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.8 0.8 0.0]; marker_styleUM(idat).m='s'; i_aerosol_processing_multi{idat}=1; append_str_timser{idat}='.txt'; idat=idat+1;
%   
%      um_str{idat}='u-an967'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='Dan, cloud scheme OFF, rain OFF';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
%   line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.8 0.0 0.0]; marker_styleUM(idat).m='s'; i_aerosol_processing_multi{idat}=1; append_str_timser{idat}='.txt'; idat=idat+1;
%   
%       um_str{idat}='u-an998'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='Annette, my switches, rain OFF';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
%   line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.0 0.0 0.8]; marker_styleUM(idat).m='o'; i_aerosol_processing_multi{idat}=1; append_str_timser{idat}='.txt'; idat=idat+1;
%    
%   
%   %%
%        um_str{idat}='u-ap417'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='Annette code, DPG switches, activation.F90 update';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
%   line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.0 0.8 0.0]; marker_styleUM(idat).m='v'; i_aerosol_processing_multi{idat}=1; append_str_timser{idat}='.txt'; idat=idat+1;
%   
%      um_str{idat}='u-ap423'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='Annette code, DPG switches, psedl=false';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
%   line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.0 0.0 0.8]; marker_styleUM(idat).m='s'; i_aerosol_processing_multi{idat}=1; append_str_timser{idat}='.txt'; idat=idat+1;
% 
%         um_str{idat}='u-ap425'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='Annette code, DPG switches, max step=10';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
%   line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.8 0.0 0.8]; marker_styleUM(idat).m='o'; i_aerosol_processing_multi{idat}=1; append_str_timser{idat}='.txt'; idat=idat+1;
%   
%      um_str{idat}='u-ap428'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='Annette code, DPG switches, Aitken ON';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
%   line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.0 0.8 0.0]; marker_styleUM(idat).m='^'; i_aerosol_processing_multi{idat}=1; append_str_timser{idat}='.txt'; idat=idat+1;
% 
%           um_str{idat}='u-ap456'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='Annette code, DPG switches,1 moment sed';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
%   line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.4 0.0 0.0]; marker_styleUM(idat).m='o'; i_aerosol_processing_multi{idat}=1; append_str_timser{idat}='.txt'; idat=idat+1;
%

%     um_str{idat}='u-ap465'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='Annette code + switches, bypass which mode';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
%  line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.0 0.0 0.4]; marker_styleUM(idat).m='^'; i_aerosol_processing_multi{idat}=1; append_str_timser{idat}='.txt'; idat=idat+1;

%  um_str{idat}='u-aq017'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='No rain sed';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
%  line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.0 0.8 0.0]; marker_styleUM(idat).m='^'; i_aerosol_processing_multi{idat}=1; append_str_timser{idat}='.txt'; idat=idat+1;

%   um_str{idat}='u-aq018'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='Per m3 No rain sed, no coarse';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
%   line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.0 0.4 0.0]; marker_styleUM(idat).m='v'; i_aerosol_processing_multi{idat}=1; append_str_timser{idat}='.txt'; idat=idat+1;

%   um_str{idat}='u-aq027'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='No rain sed, no coarse, no pos';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
%   line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.0 0.0 0.8]; marker_styleUM(idat).m='o'; i_aerosol_processing_multi{idat}=1; append_str_timser{idat}='.txt'; idat=idat+1;
% 
%   um_str{idat}='u-aq028'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='No rain sed, no coarse, no pos, no tidy';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
%   line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.0 0.0 0.4]; marker_styleUM(idat).m='s'; i_aerosol_processing_multi{idat}=1; append_str_timser{idat}='.txt'; idat=idat+1;
% 
%   um_str{idat}='u-aq034'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='No rain sed, no coarse, no pos, no tidy, no asedl';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
%   line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.0 0.0 0.4]; marker_styleUM(idat).m='d'; i_aerosol_processing_multi{idat}=1; append_str_timser{idat}='.txt'; idat=idat+1;
% 
%   um_str{idat}='u-aq036'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='No rain sed, no coarse, no asedl';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
%   line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.8 0.0 0.0]; marker_styleUM(idat).m='^'; i_aerosol_processing_multi{idat}=1; append_str_timser{idat}='.txt'; idat=idat+1;

%   um_str{idat}='u-aq164'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='Per m3 Coarse mode on, aerosol present';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
%   line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.0 0.0 0.4]; marker_styleUM(idat).m='s'; i_aerosol_processing_multi{idat}=1; append_str_timser{idat}='.txt'; idat=idat+1;

%  um_str{idat}='u-aq165'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='Coarse, Aitken on, aerosol in coarse only';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
%  line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.0 0.0 0.4]; marker_styleUM(idat).m='d'; i_aerosol_processing_multi{idat}=1; append_str_timser{idat}='.txt'; idat=idat+1;

%  um_str{idat}='u-aq167'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='Coarse, Aitken on, aerosol in both';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
%  line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.8 0.0 0.0]; marker_styleUM(idat).m='^'; i_aerosol_processing_multi{idat}=1; append_str_timser{idat}='.txt'; idat=idat+1;

%um_str{idat}='u-aq221'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='Passive, accum only';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
%line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.0 0.8 0.0]; marker_styleUM(idat).m='^'; i_aerosol_processing_multi{idat}=1; append_str_timser{idat}='.txt'; idat=idat+1;

%um_str{idat}='u-aq251'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='l-mr-physics1=false, passive, accum only';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
%line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.4 0.0 0.0]; marker_styleUM(idat).m='o'; i_aerosol_processing_multi{idat}=1; append_str_timser{idat}='.txt'; idat=idat+1;

%um_str{idat}='u-aq254'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='Per kg advection, psedl ON';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
%line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.0 0.0 0.4]; marker_styleUM(idat).m='^'; i_aerosol_processing_multi{idat}=1; append_str_timser{idat}='.txt'; idat=idat+1;

% um_str{idat}='u-aq483'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='Per kg advection, psedl OFF, coarse mode';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
% line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.0 0.0 0.4]; marker_styleUM(idat).m='^'; i_aerosol_processing_multi{idat}=1; append_str_timser{idat}='.txt'; idat=idat+1;

% um_str{idat}='u-aq485'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='Per m3 advection, psedl OFF, coarse mode';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
% line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.4 0.0 0.0]; marker_styleUM(idat).m='o'; i_aerosol_processing_multi{idat}=1; append_str_timser{idat}='.txt'; idat=idat+1;

 um_str{idat}='u-aq608'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='Per kg, 50% act, psedl OFF, coarse mode';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
 line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.0 0.8 0.0]; marker_styleUM(idat).m='^'; i_aerosol_processing_multi{idat}=1; append_str_timser{idat}='.txt'; idat=idat+1;
% 
um_str{idat}='u-aq709'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='Per kg, 50% act, no sec act, psedl OFF, coarse mode';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.0 0.0 0.8]; marker_styleUM(idat).m='v'; i_aerosol_processing_multi{idat}=1; append_str_timser{idat}='.txt'; idat=idat+1;

% um_str{idat}='u-aq922'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='Bug fix, per kg, ARG, coarse mode';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
% line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.8 0.0 0.0]; marker_styleUM(idat).m='v'; i_aerosol_processing_multi{idat}=1; append_str_timser{idat}='.txt'; idat=idat+1;
% 
% um_str{idat}='u-aq924'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='Bug fix, per kg, ARG, coarse OFF';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
% line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.0 0.8 0.0]; marker_styleUM(idat).m='v'; i_aerosol_processing_multi{idat}=1; append_str_timser{idat}='.txt'; idat=idat+1;
% 
% um_str{idat}='u-aq973'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='Bug fix, per kg, 50% act, no sec, coarse OFF';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
%  line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.4 0.0 0.0]; marker_styleUM(idat).m='o'; i_aerosol_processing_multi{idat}=1; append_str_timser{idat}='.txt'; idat=idat+1;
% 
% um_str{idat}='u-aq974'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='Bug fix, per kg, 98% act, no sec, coarse mode ON';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
% line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.0 0.0 0.4]; marker_styleUM(idat).m='s'; i_aerosol_processing_multi{idat}=1; append_str_timser{idat}='.txt'; idat=idat+1;
%  
um_str{idat}='u-aq993'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='Bug fix,per kg,50% act,no sec,coarse mode ON,no aero';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.0 0.4 0.0]; marker_styleUM(idat).m='v'; i_aerosol_processing_multi{idat}=1; append_str_timser{idat}='.txt'; idat=idat+1;

% um_str{idat}='u-aq994'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='Bug fix, per kg, ARG act, no sec, coarse mode ON';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
% line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.0 0.0 0.4]; marker_styleUM(idat).m='o'; i_aerosol_processing_multi{idat}=1; append_str_timser{idat}='.txt'; idat=idat+1;



%um_str{idat}='u-ar023'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='Bug fix OFF, per kg, 50% act, no sec, coarse OFF, no aero';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
%line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.0 0.0 0.4]; marker_styleUM(idat).m='o'; i_aerosol_processing_multi{idat}=1; append_str_timser{idat}='.txt'; idat=idat+1;

% um_str{idat}='u-ar021'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='Bug fix,per kg,50% act,no sec,coarse OFF,no aero';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
% line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.4 0.4 0.0]; marker_styleUM(idat).m='s'; i_aerosol_processing_multi{idat}=1; append_str_timser{idat}='.txt'; idat=idat+1;

%um_str{idat}='u-ar024'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='Bug fix, per kg, 50% act, sec ON, coarse OFF, no aero';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
%line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.0 0.0 0.8]; marker_styleUM(idat).m='v'; i_aerosol_processing_multi{idat}=1; append_str_timser{idat}='.txt'; idat=idat+1;

%um_str{idat}='u-ar033'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='Bug fix + dm=0, per kg, ARG, coarse mode';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
% line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.4 0.0 0.0]; marker_styleUM(idat).m='o'; i_aerosol_processing_multi{idat}=1; append_str_timser{idat}='.txt'; idat=idat+1;

%um_str{idat}='u-ar036'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='Bug fix,per kg,ARG,coarse mode,no aero,evap to coarse';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
%line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.0 0.0 0.4]; marker_styleUM(idat).m='s'; i_aerosol_processing_multi{idat}=1; append_str_timser{idat}='.txt'; idat=idat+1;
%Actually used the namelist with coarse aerosol - some aerosol is
%transferred to coarse through evap

%um_str{idat}='u-ar086'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='Bug fix, per kg, ARG, coarse mode as accum';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
%line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.0 0.0 0.4]; marker_styleUM(idat).m='s'; i_aerosol_processing_multi{idat}=1; append_str_timser{idat}='.txt'; idat=idat+1;

%  um_str{idat}='u-ar087'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='Bug fix OFF,per kg,50% act,no sec, coarse mode ON,no aero';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
%  line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.8 0.0 0.0]; marker_styleUM(idat).m='v'; i_aerosol_processing_multi{idat}=1; append_str_timser{idat}='.txt'; idat=idat+1;
% % 
  um_str{idat}='u-ar089'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='Bug fix,per kg,50% act,sec ON,coarse mode ON,no aero';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
  line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.0 0.8 0.0]; marker_styleUM(idat).m='v'; i_aerosol_processing_multi{idat}=1; append_str_timser{idat}='.txt'; idat=idat+1;
% % 
% um_str{idat}='u-ar056'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='Bug fix + dm=0, per kg, ARG, coarse mode, tapered nml';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
% line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.0 0.0 0.4]; marker_styleUM(idat).m='s'; i_aerosol_processing_multi{idat}=1; append_str_timser{idat}='.txt'; idat=idat+1;

% --- u-ar268 useful as a reference for mass changes with no activation ---
%  um_str{idat}='u-ar268'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='NO ACT, Bug fix,per kg,50% act, coarse mode act';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
%  line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.8 0.0 0.0]; marker_styleUM(idat).m='^'; i_aerosol_processing_multi{idat}=1; append_str_timser{idat}='.txt'; idat=idat+1;

% --- u-ar525 useful as a reference for mass changes with no activation ---
  um_str{idat}='u-ar525'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='NO ACT, Bug fix,per kg,50% act, coarse mode act';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
  line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.8 0.0 0.0]; marker_styleUM(idat).m='^'; i_aerosol_processing_multi{idat}=1; append_str_timser{idat}='.txt'; idat=idat+1;
  
% --- u-ar525 except with coarse mode set to zero (manually in .txt files, copy of aitken) ---
  um_str{idat}='u-ar525-zero-coarse'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='NO ACT, Bug fix,per kg,50% act, coarse mode act';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
  line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.0 0.0 0.8]; marker_styleUM(idat).m='^'; i_aerosol_processing_multi{idat}=1; append_str_timser{idat}='.txt'; idat=idat+1;
    
  
%  um_str{idat}='u-ar269'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='Bug fix,per kg,50% act,no sec,accum mode act';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
%  line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.0 0.0 0.8]; marker_styleUM(idat).m='^'; i_aerosol_processing_multi{idat}=1; append_str_timser{idat}='.txt'; idat=idat+1;

%  um_str{idat}='u-ar364'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='u-ar364';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
%  line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.0 0.4 0.0]; marker_styleUM(idat).m='v'; i_aerosol_processing_multi{idat}=1; append_str_timser{idat}='.txt'; idat=idat+1;

%  um_str{idat}='u-ar365'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='u-ar365';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
%  line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.0 0.0 0.4]; marker_styleUM(idat).m='o'; i_aerosol_processing_multi{idat}=1; append_str_timser{idat}='.txt'; idat=idat+1;

% um_str{idat}='u-ar557'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='Per kg nml, no act';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
% line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.0 0.0 0.8]; marker_styleUM(idat).m='o'; i_aerosol_processing_multi{idat}=1; append_str_timser{idat}='.txt'; idat=idat+1;

%  um_str{idat}='u-as088'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='Per kg nml, act of coarse only';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
%  line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.0 0.4 0.0]; marker_styleUM(idat).m='v'; i_aerosol_processing_multi{idat}=1; append_str_timser{idat}='.txt'; idat=idat+1;
% 
%  um_str{idat}='u-as093'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='Per kg nml, act of accum only';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
%  line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.0 0.0 0.4]; marker_styleUM(idat).m='o'; i_aerosol_processing_multi{idat}=1; append_str_timser{idat}='.txt'; idat=idat+1;
% 
%  um_str{idat}='u-as420'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='Per kg nml, NO ACT';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
%  line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.8 0.0 0.0]; marker_styleUM(idat).m='o'; i_aerosol_processing_multi{idat}=1; append_str_timser{idat}='.txt'; idat=idat+1;
% 
%  um_str{idat}='u-as421'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='Per kg nml, NO ACT %=0';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
%  line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.0 0.0 0.8]; marker_styleUM(idat).m='o'; i_aerosol_processing_multi{idat}=1; append_str_timser{idat}='.txt'; idat=idat+1;


            case 'per kg namelist'
 
        
%% Old per kg namlelist runs   :-

%  um_str{idat}='u-aq608'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='Per kg, 50% act, psedl OFF, coarse mode';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
%  line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.0 0.8 0.0]; marker_styleUM(idat).m='^'; i_aerosol_processing_multi{idat}=1; append_str_timser{idat}='.txt'; idat=idat+1;
% % 
% um_str{idat}='u-aq709'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='Per kg, 50% act, no sec act, psedl OFF, coarse mode';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
% line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.0 0.0 0.8]; marker_styleUM(idat).m='v'; i_aerosol_processing_multi{idat}=1; append_str_timser{idat}='.txt'; idat=idat+1;
% 
% um_str{idat}='u-aq993'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='Bug fix,per kg,50% act,no sec,coarse mode ON,no aero';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
% line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.0 0.4 0.0]; marker_styleUM(idat).m='v'; i_aerosol_processing_multi{idat}=1; append_str_timser{idat}='.txt'; idat=idat+1;
% 
%   um_str{idat}='u-ar089'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='Bug fix,per kg,50% act,sec ON,coarse mode ON,no aero';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
%   line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.0 0.8 0.0]; marker_styleUM(idat).m='v'; i_aerosol_processing_multi{idat}=1; append_str_timser{idat}='.txt'; idat=idat+1;
% % % 

% --- u-ar525 useful as a reference for mass changes with no activation ---
%  um_str{idat}='u-ar525'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='NO ACT, Bug fix,per kg,50% act, coarse mode act';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
%  line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.8 0.0 0.0]; marker_styleUM(idat).m='^'; i_aerosol_processing_multi{idat}=1; append_str_timser{idat}='.txt'; idat=idat+1;   

%% New per kg namlelist runs   :-
    
%  um_str{idat}='u-as088'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='Coarse to coarse';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
%  line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.0 0.4 0.0]; marker_styleUM(idat).m='v'; i_aerosol_processing_multi{idat}=1; append_str_timser{idat}='.txt'; idat=idat+1;
% 
%  um_str{idat}='u-as093'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='Accum to accum';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
%  line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.0 0.0 0.4]; marker_styleUM(idat).m='o'; i_aerosol_processing_multi{idat}=1; append_str_timser{idat}='.txt'; idat=idat+1;
% 
% % um_str{idat}='u-as420'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='No aero in nml';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
% % line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.8 0.0 0.0]; marker_styleUM(idat).m='o'; i_aerosol_processing_multi{idat}=1; append_str_timser{idat}='.txt'; idat=idat+1;
% 
%  um_str{idat}='u-as421'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='No activation';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
%  line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.0 0.0 0.8]; marker_styleUM(idat).m='o'; i_aerosol_processing_multi{idat}=1; append_str_timser{idat}='.txt'; idat=idat+1;

 um_str{idat}='u-at545'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='Act both, return to accum';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
 line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.4 0.4 0.0]; marker_styleUM(idat).m='o'; i_aerosol_processing_multi{idat}=1; append_str_timser{idat}='.txt'; idat=idat+1;

 um_str{idat}='u-at546'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='Act both, return to coarse';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
 line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.0 0.4 0.4]; marker_styleUM(idat).m='o'; i_aerosol_processing_multi{idat}=1; append_str_timser{idat}='.txt'; idat=idat+1;

% um_str{idat}='u-av650'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='Act accum, return to accum, zero coarse';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
% line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.4 0.0 0.0]; marker_styleUM(idat).m='o'; i_aerosol_processing_multi{idat}=1; append_str_timser{idat}='.txt'; idat=idat+1;

% um_str{idat}='u-av652'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='No act, zero coarse';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
% line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.8 0.8 0.0]; marker_styleUM(idat).m='o'; i_aerosol_processing_multi{idat}=1; append_str_timser{idat}='.txt'; idat=idat+1;

 um_str{idat}='u-av714'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='Act coarse, return to accum';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
 line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.4 0.0 0.0]; marker_styleUM(idat).m='o'; i_aerosol_processing_multi{idat}=1; append_str_timser{idat}='.txt'; idat=idat+1;

 um_str{idat}='u-av715'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='Act accum, return to coarse';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
 line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.8 0.8 0.0]; marker_styleUM(idat).m='o'; i_aerosol_processing_multi{idat}=1; append_str_timser{idat}='.txt'; idat=idat+1;
 
 

 um_str{idat}='u-av955'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='Act both, NO return';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
 line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.0 0.4 0.4]; marker_styleUM(idat).m='s'; i_aerosol_processing_multi{idat}=1; append_str_timser{idat}='.txt'; idat=idat+1;

  um_str{idat}='u-av953'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='Act coarse, NO return';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
 line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.4 0.0 0.0]; marker_styleUM(idat).m='s'; i_aerosol_processing_multi{idat}=1; append_str_timser{idat}='.txt'; idat=idat+1;

 um_str{idat}='u-av952'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='Act accum, NO return';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
 line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.8 0.8 0.0]; marker_styleUM(idat).m='s'; i_aerosol_processing_multi{idat}=1; append_str_timser{idat}='.txt'; idat=idat+1;

 um_str{idat}='u-av967'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='Act both, return to accum, tracer for act';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
 line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.4 0.4 0.0]; marker_styleUM(idat).m='o'; i_aerosol_processing_multi{idat}=1; append_str_timser{idat}='.txt'; idat=idat+1;

  um_str{idat}='u-aw107'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='Act both, return to accum, dmass test';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
 line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.4 0.4 0.0]; marker_styleUM(idat).m='s'; i_aerosol_processing_multi{idat}=1; append_str_timser{idat}='.txt'; idat=idat+1;

 um_str{idat}='u-aw110'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='Act both, NO return, NO tidy';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
 line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.8 0.8 0.0]; marker_styleUM(idat).m='s'; i_aerosol_processing_multi{idat}=1; append_str_timser{idat}='.txt'; idat=idat+1;

        end

    case '12th Nov case, cloud scheme test v10.8'
        
        dirUM2{1}='/home/disk/eos8/d.grosvenor/UM/12Nov2008_Boutle/';
        dirUM2{2}='/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/'; 
        
        idat=1;
        
        um_str{idat}='u-av386'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='CF=1 for mphys';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
            line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.0 0.4 0.4]; marker_styleUM(idat).m='o'; i_aerosol_processing_multi{idat}=1; append_str_timser{idat}='.txt'; idat=idat+1;
        um_str{idat}='u-av387'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='CF=1 for rad';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
         line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.4 0.0 0.4]; marker_styleUM(idat).m='v'; i_aerosol_processing_multi{idat}=1; append_str_timser{idat}='.txt'; idat=idat+1;
        um_str{idat}='u-av399'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='Full cloud scheme';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
            line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.4 0.4 0.0]; marker_styleUM(idat).m='s'; i_aerosol_processing_multi{idat}=1; append_str_timser{idat}='.txt'; idat=idat+1;
        um_str{idat}='u-av410'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='Cloud scheme OFF';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
            line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.8 0.0 0.0]; marker_styleUM(idat).m='^'; i_aerosol_processing_multi{idat}=1; append_str_timser{idat}='.txt'; idat=idat+1;


    case '12th Nov case, as of May 2016 rain OFF multi-dirUM'
        
        dirUM2{1}='/home/disk/eos8/d.grosvenor/UM/12Nov2008_Boutle/';                
        dirUM2{2}='/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/';        

        idat=1;
        
%Rain off runs        
       um_str{idat}='u-ai120'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='CASIMv10.3-Ndvar-no-auto';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
       line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.8 0.4 0]; marker_styleUM(idat).m='^'; idat=idat+1;
        
       um_str{idat}='u-ai117'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='CASIMv10.3-Ndvar-no-auto-proc';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
       line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.8 0.4 0]; marker_styleUM(idat).m='^'; idat=idat+1;
  
       um_str{idat}='u-ai121'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='CASIMv10.3-Ndvar-no-auto-proc-cloudSchemeOFF';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
       line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.4 0.4 0.4]; marker_styleUM(idat).m='^'; idat=idat+1;
       
         um_str{idat}='u-ai123'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='CASIMv10.3-Ndvar-no-auto-cloudSchemeOFF';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
         line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.4 0.4 0.4]; marker_styleUM(idat).m='^'; idat=idat+1;             
% 
        um_str{idat}='u-ai397'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='CASIMv10.3-Ndvar-no-auto-proc-cloudSchemeOFF-surf-source';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
        line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0 0.4 0.8]; marker_styleUM(idat).m='^'; idat=idat+1;        

        um_str{idat}='u-ai417'; fileUM{idat} = [um_str{idat} '/' um_str{idat} '_VAR_NAME_VOCALS_1p0_L70_ukv_.pp.nc']; dirUM{idat}=dirUM2{2}; labs_UM(idat).l ='CASIMv10.3-Ndvar-Annette';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
        line_patternUM(idat).p= '-';  line_colourUM(idat).c=[0.4 0.8 0.4]; marker_styleUM(idat).m='^'; idat=idat+1;



        
    case '12th Nov case, simple for Steve Abel'
        %% Newest runs May 2016


        dirUM='/home/disk/eos8/d.grosvenor/UM/26thOct_POC/';
        dirUM='/home/disk/eos8/d.grosvenor/UM/12Nov2008_Boutle/';

        idat=1;
        fileUM{idat} = '/xmmz-u/xmmzu_VAR_NAME_.pp.nc'; labs_UM(idat).l = 'CASIM-Ndvar';  flag{idat} = 'load_mat'; fileUM_rho{idat} = '/xmmz-u/xmmzu_rho_.pp.nc';pole_lat=70; pole_lon=284;
        line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.4 0.4 0.4]; marker_styleUM(idat).m='d'; idat=idat+1;
%        fileUM{idat} = '/xmmz-v/xmmzv_VAR_NAME_.pp.nc'; labs_UM(idat).l ='CASIM-Ndvar-0.1';  flag{idat} = 'load_mat'; fileUM_rho{idat} = '/xmmz-v/xmmzv_rho_.pp.nc'; pole_lat=70; pole_lon=284;
%        line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0 0 0]; marker_styleUM(idat).m='v'; idat=idat+1;
%        fileUM{idat} = '/xmmz-w/xmmzw_VAR_NAME_.pp.nc'; labs_UM(idat).l ='CASIM-Ndvar-10';  flag{idat} = 'load_mat'; fileUM_rho{idat} = '/xmmz-v/xmmzv_rho_.pp.nc'; pole_lat=70; pole_lon=284;
%        line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.8 0.4 0]; marker_styleUM(idat).m='^'; idat=idat+1;
        %fileUM{idat} = '/xlhg-v/xlhgv_VAR_NAME_.pp.nc'; labs_UM(idat).l = 'CASIM-Ndvar-0.025';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
        %        line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0 0.8 0]; marker_styleUM(idat).m='v'; idat=idat+1;
        fileUM{idat} = '/xmmz-x/xmmzx_VAR_NAME_.pp.nc'; labs_UM(idat).l = 'CASIM-Ndvar-0.025';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
        line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0 0.8 0]; marker_styleUM(idat).m='o'; idat=idat+1;
%        fileUM{idat} = '/xmmz-n/xmmzn_VAR_NAME_.pp.nc'; labs_UM(idat).l = 'Old-mphys';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
 %       line_patternUM(idat).p= '--';  line_colourUM(idat).c=[1 0 0]; marker_styleUM(idat).m='s'; idat=idat+1;
        %fileUM{idat} = '/xmmz-l/xmmzl_VAR_NAME_.pp.nc'; labs_UM(idat).l ='CASIM-Ndvar-RHcrit0.999';  flag{idat} = 'load_mat'; fileUM_rho{idat} = '/xmmz-v/xmmzv_rho_.pp.nc'; pole_lat=70; pole_lon=284;
        %    line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0 0 0];
        %    marker_styleUM(idat).m='v'; idat=idat+1;
    
    case '12th Nov case ACI, as of May 2016'
        %% Newest runs May 2016


        dirUM='/home/disk/eos8/d.grosvenor/UM/26thOct_POC/';
        dirUM='/home/disk/eos8/d.grosvenor/UM/12Nov2008_Boutle/';

        idat=1;
        fileUM{idat} = '/xmmz-u/xmmzu_VAR_NAME_.pp.nc'; labs_UM(idat).l = 'CASIM-Ndvar';  flag{idat} = 'load_mat'; fileUM_rho{idat} = '/xmmz-u/xmmzu_rho_.pp.nc';pole_lat=70; pole_lon=284;
        line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.4 0.4 0.4]; marker_styleUM(idat).m='d'; idat=idat+1;
        fileUM{idat} = '/xmmz-v/xmmzv_VAR_NAME_.pp.nc'; labs_UM(idat).l ='CASIM-Ndvar-0.1';  flag{idat} = 'load_mat'; fileUM_rho{idat} = '/xmmz-v/xmmzv_rho_.pp.nc'; pole_lat=70; pole_lon=284;
        line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0 0 0]; marker_styleUM(idat).m='v'; idat=idat+1;
        fileUM{idat} = '/xmmz-w/xmmzw_VAR_NAME_.pp.nc'; labs_UM(idat).l ='CASIM-Ndvar-10';  flag{idat} = 'load_mat'; fileUM_rho{idat} = '/xmmz-v/xmmzv_rho_.pp.nc'; pole_lat=70; pole_lon=284;
        line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.8 0.4 0]; marker_styleUM(idat).m='^'; idat=idat+1;
        %fileUM{idat} = '/xlhg-v/xlhgv_VAR_NAME_.pp.nc'; labs_UM(idat).l = 'CASIM-Ndvar-0.025';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
        %        line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0 0.8 0]; marker_styleUM(idat).m='v'; idat=idat+1;
        fileUM{idat} = '/xmmz-x/xmmzx_VAR_NAME_.pp.nc'; labs_UM(idat).l = 'CASIM-Ndvar-0.025';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
        line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0 0.8 0]; marker_styleUM(idat).m='o'; idat=idat+1;
%        fileUM{idat} = '/xmmz-n/xmmzn_VAR_NAME_.pp.nc'; labs_UM(idat).l = 'Old-mphys';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
%        line_patternUM(idat).p= '--';  line_colourUM(idat).c=[1 0 0]; marker_styleUM(idat).m='s'; idat=idat+1;
        %fileUM{idat} = '/xmmz-l/xmmzl_VAR_NAME_.pp.nc'; labs_UM(idat).l ='CASIM-Ndvar-RHcrit0.999';  flag{idat} = 'load_mat'; fileUM_rho{idat} = '/xmmz-v/xmmzv_rho_.pp.nc'; pole_lat=70; pole_lon=284;
        %    line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0 0 0]; marker_styleUM(idat).m='v'; idat=idat+1;

 case '12th Nov case ACI, as of May 2016 with sed runs'
        %% Newest runs May 2016


        dirUM='/home/disk/eos8/d.grosvenor/UM/26thOct_POC/';
        dirUM='/home/disk/eos8/d.grosvenor/UM/12Nov2008_Boutle/';

        idat=1;
        fileUM{idat} = '/xmmz-u/xmmzu_VAR_NAME_.pp.nc'; labs_UM(idat).l = 'CASIM-Ndvar';  flag{idat} = 'load_mat'; fileUM_rho{idat} = '/xmmz-u/xmmzu_rho_.pp.nc';pole_lat=70; pole_lon=284;
        line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.4 0.4 0.4]; marker_styleUM(idat).m='d'; idat=idat+1;
        fileUM{idat} = '/xmmz-v/xmmzv_VAR_NAME_.pp.nc'; labs_UM(idat).l ='CASIM-Ndvar-0.1';  flag{idat} = 'load_mat'; fileUM_rho{idat} = '/xmmz-v/xmmzv_rho_.pp.nc'; pole_lat=70; pole_lon=284;
        line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0 0 0]; marker_styleUM(idat).m='v'; idat=idat+1;
        fileUM{idat} = '/xmmz-w/xmmzw_VAR_NAME_.pp.nc'; labs_UM(idat).l ='CASIM-Ndvar-10';  flag{idat} = 'load_mat'; fileUM_rho{idat} = '/xmmz-v/xmmzv_rho_.pp.nc'; pole_lat=70; pole_lon=284;
        line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.8 0.4 0]; marker_styleUM(idat).m='^'; idat=idat+1;
        %fileUM{idat} = '/xlhg-v/xlhgv_VAR_NAME_.pp.nc'; labs_UM(idat).l = 'CASIM-Ndvar-0.025';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
        %        line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0 0.8 0]; marker_styleUM(idat).m='v'; idat=idat+1;
        fileUM{idat} = '/xmmz-x/xmmzx_VAR_NAME_.pp.nc'; labs_UM(idat).l = 'CASIM-Ndvar-0.025';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
        line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0 0.8 0]; marker_styleUM(idat).m='o'; idat=idat+1;
%        fileUM{idat} = '/xmmz-n/xmmzn_VAR_NAME_.pp.nc'; labs_UM(idat).l = 'Old-mphys';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
%        line_patternUM(idat).p= '--';  line_colourUM(idat).c=[1 0 0]; marker_styleUM(idat).m='s'; idat=idat+1;
        %fileUM{idat} = '/xmmz-l/xmmzl_VAR_NAME_.pp.nc'; labs_UM(idat).l ='CASIM-Ndvar-RHcrit0.999';  flag{idat} = 'load_mat'; fileUM_rho{idat} = '/xmmz-v/xmmzv_rho_.pp.nc'; pole_lat=70; pole_lon=284;
        %    line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0 0 0]; marker_styleUM(idat).m='v'; idat=idat+1;
        fileUM{idat} = '/xmmz-p/xmmzp_VAR_NAME_.pp.nc'; labs_UM(idat).l = 'CASIM-Ndvar-sed';  flag{idat} = 'load_mat'; fileUM_rho{idat} = '';pole_lat=70; pole_lon=284;
        line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.4 0.4 0.4]; marker_styleUM(idat).m='s'; idat=idat+1;    
        fileUM{idat} = '/xmmz-q/xmmzq_VAR_NAME_.pp.nc'; labs_UM(idat).l = 'CASIM-Ndvar-10-sed';  flag{idat} = 'load_mat'; fileUM_rho{idat} = '';pole_lat=70; pole_lon=284;
        line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.8 0.4 0]; marker_styleUM(idat).m='s'; idat=idat+1;                   

    case 'Iceland_9day_runs_Nov2016'

        dirUM='/home/disk/eos10/d.grosvenor/UM/Iceland/';        

        idat=1;
       fileUM{idat} = '/u-af148/u-af148_VAR_NAME_Iceland_4p0_L70_ukv_.pp.nc'; labs_UM(idat).l = 'Low aerosol, volcano OFF';  flag{idat} = 'load_mat'; pole_lat=25; pole_lon=165;
       line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.4 0.4 0.4]; marker_styleUM(idat).m='d'; idat=idat+1;
       fileUM{idat} = '/u-af178/u-af178_VAR_NAME_Iceland_4p0_L70_ukv_.pp.nc'; labs_UM(idat).l ='High aerosol, volcano ON';  flag{idat} = 'load_mat'; pole_lat=25; pole_lon=165;
       line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0 0 0]; marker_styleUM(idat).m='v'; idat=idat+1;
       fileUM{idat} = '/u-af249/u-af249_VAR_NAME_Iceland_4p0_L70_ukv_.pp.nc'; labs_UM(idat).l ='Low aerosol, volcano ON';  flag{idat} = 'load_mat';  pole_lat=25; pole_lon=165;
       line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.8 0.4 0]; marker_styleUM(idat).m='^'; idat=idat+1;
        fileUM{idat} = '/u-af250/u-af250_VAR_NAME_Iceland_4p0_L70_ukv_.pp.nc'; labs_UM(idat).l = 'High aerosol, volcano OFF';  flag{idat} = 'load_mat'; pole_lat=25; pole_lon=165;
        line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0 0.8 0]; marker_styleUM(idat).m='v'; idat=idat+1;

        
    case 'Iceland_9day_runs_Nov2016 adhoc'

        dirUM='/home/disk/eos10/d.grosvenor/UM/Iceland/';        

        idat=1;
%       fileUM{idat} = '/u-af148/u-af148_VAR_NAME_Iceland_4p0_L70_ukv_.pp.nc'; labs_UM(idat).l = 'Low aerosol, volcano OFF';  flag{idat} = 'load_mat'; pole_lat=25; pole_lon=165;
%       line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.4 0.4 0.4]; marker_styleUM(idat).m='d'; idat=idat+1;
%       fileUM{idat} = '/u-af178/u-af178_VAR_NAME_Iceland_4p0_L70_ukv_.pp.nc'; labs_UM(idat).l ='High aerosol, volcano ON';  flag{idat} = 'load_mat'; pole_lat=25; pole_lon=165;
%       line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0 0 0]; marker_styleUM(idat).m='v'; idat=idat+1;
       fileUM{idat} = '/u-af249/u-af249_VAR_NAME_Iceland_4p0_L70_ukv_.pp.nc'; labs_UM(idat).l ='Low aerosol, volcano ON';  flag{idat} = 'load_mat';  pole_lat=25; pole_lon=165;
       line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.8 0.4 0]; marker_styleUM(idat).m='^'; idat=idat+1;
%        fileUM{idat} = '/u-af250/u-af250_VAR_NAME_Iceland_4p0_L70_ukv_.pp.nc'; labs_UM(idat).l = 'High aerosol, volcano OFF';  flag{idat} = 'load_mat'; pole_lat=25; pole_lon=165;
%        line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0 0.8 0]; marker_styleUM(idat).m='v'; idat=idat+1;


     case 'Iceland_9day_runs_Nov2016 - NEW runs'

        dirUM='/home/disk/eos10/d.grosvenor/UM/Iceland/';        

        idat=1;
        fileUM{idat} = '/u-ai864/u-ai864_VAR_NAME_Iceland_4p0_L70_ukv_.pp.nc'; labs_UM(idat).l = 'Low aerosol, volcano OFF';  flag{idat} = 'load_mat'; pole_lat=25; pole_lon=165;
        line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0 0.8 0]; marker_styleUM(idat).m='v'; idat=idat+1;
        fileUM{idat} = '/u-ai866/u-ai866_VAR_NAME_Iceland_4p0_L70_ukv_.pp.nc'; labs_UM(idat).l = 'Low aerosol, volcano ON';  flag{idat} = 'load_mat'; pole_lat=25; pole_lon=165;
        line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0 0.8 0]; marker_styleUM(idat).m='v'; idat=idat+1;       
        
    case 'Iceland_9day_runs_Nov2016 - volcanoes ON only'

        dirUM='/home/disk/eos10/d.grosvenor/UM/Iceland/';        

        idat=1;
%       fileUM{idat} = '/u-af148/u-af148_VAR_NAME_Iceland_4p0_L70_ukv_.pp.nc'; labs_UM(idat).l = 'Low aerosol, volcano OFF';  flag{idat} = 'load_mat'; pole_lat=25; pole_lon=165;
%       line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.4 0.4 0.4]; marker_styleUM(idat).m='d'; idat=idat+1;
       fileUM{idat} = '/u-af178/u-af178_VAR_NAME_Iceland_4p0_L70_ukv_.pp.nc'; labs_UM(idat).l ='High aerosol, volcano ON';  flag{idat} = 'load_mat'; pole_lat=25; pole_lon=165;
       line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0 0 0]; marker_styleUM(idat).m='v'; idat=idat+1;
       fileUM{idat} = '/u-af249/u-af249_VAR_NAME_Iceland_4p0_L70_ukv_.pp.nc'; labs_UM(idat).l ='Low aerosol, volcano ON';  flag{idat} = 'load_mat';  pole_lat=25; pole_lon=165;
       line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.8 0.4 0]; marker_styleUM(idat).m='^'; idat=idat+1;
%        fileUM{idat} = '/u-af250/u-af250_VAR_NAME_Iceland_4p0_L70_ukv_.pp.nc'; labs_UM(idat).l = 'High aerosol, volcano OFF';  flag{idat} = 'load_mat'; pole_lat=25; pole_lon=165;
%        line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0 0.8 0]; marker_styleUM(idat).m='v'; idat=idat+1;
       

    case 'Iceland_9day_runs_Nov2016_low_background_only'

        dirUM='/home/disk/eos10/d.grosvenor/UM/Iceland/';        

        idat=1;
       fileUM{idat} = '/u-af148/u-af148_VAR_NAME_Iceland_4p0_L70_ukv_.pp.nc'; labs_UM(idat).l = 'Low aerosol, volcano OFF';  flag{idat} = 'load_mat'; pole_lat=25; pole_lon=165;
       line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.4 0.4 0.4]; marker_styleUM(idat).m='d'; idat=idat+1;
%       fileUM{idat} = '/u-af178/u-af178_VAR_NAME_Iceland_4p0_L70_ukv_.pp.nc'; labs_UM(idat).l ='High aerosol, volcano ON';  flag{idat} = 'load_mat'; pole_lat=25; pole_lon=165;
%       line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0 0 0]; marker_styleUM(idat).m='v'; idat=idat+1;
       fileUM{idat} = '/u-af249/u-af249_VAR_NAME_Iceland_4p0_L70_ukv_.pp.nc'; labs_UM(idat).l ='Low aerosol, volcano ON';  flag{idat} = 'load_mat';  pole_lat=25; pole_lon=165;
       line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.8 0.4 0]; marker_styleUM(idat).m='^'; idat=idat+1;
%        fileUM{idat} = '/u-af250/u-af250_VAR_NAME_Iceland_4p0_L70_ukv_.pp.nc'; labs_UM(idat).l = 'High aerosol, volcano OFF';  flag{idat} = 'load_mat'; pole_lat=25; pole_lon=165;
%        line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0 0.8 0]; marker_styleUM(idat).m='v'; idat=idat+1;

case 'Iceland_9day_runs_Nov2016_low_background_volcano_ON_only'

        dirUM='/home/disk/eos10/d.grosvenor/UM/Iceland/';        

        idat=1;
%       fileUM{idat} = '/u-af148/u-af148_VAR_NAME_Iceland_4p0_L70_ukv_.pp.nc'; labs_UM(idat).l = 'Low aerosol, volcano OFF';  flag{idat} = 'load_mat'; pole_lat=25; pole_lon=165;
%       line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.4 0.4 0.4]; marker_styleUM(idat).m='d'; idat=idat+1;
%       fileUM{idat} = '/u-af178/u-af178_VAR_NAME_Iceland_4p0_L70_ukv_.pp.nc'; labs_UM(idat).l ='High aerosol, volcano ON';  flag{idat} = 'load_mat'; pole_lat=25; pole_lon=165;
%       line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0 0 0]; marker_styleUM(idat).m='v'; idat=idat+1;
       fileUM{idat} = '/u-af249/u-af249_VAR_NAME_Iceland_4p0_L70_ukv_.pp.nc'; labs_UM(idat).l ='Low aerosol, volcano ON';  flag{idat} = 'load_mat';  pole_lat=25; pole_lon=165;
       line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.8 0.4 0]; marker_styleUM(idat).m='^'; idat=idat+1;
%        fileUM{idat} = '/u-af250/u-af250_VAR_NAME_Iceland_4p0_L70_ukv_.pp.nc'; labs_UM(idat).l = 'High aerosol, volcano OFF';  flag{idat} = 'load_mat'; pole_lat=25; pole_lon=165;
%        line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0 0.8 0]; marker_styleUM(idat).m='v'; idat=idat+1;


    case 'UKESM'
        %Copied these from calc_LWP_multi - might need altering to work
        %UKESM eval May 2016
        % fileUM{idat} = '/u-ab642/u-ab642_Nd_times_LYR_.pp.nc'; labs_UM(idat).l = 'u-ab642';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284; idat=idat+1;
        % fileUM{idat} = '/u-ac043/u-ac043_Nd_times_LYR_.pp.nc'; labs_UM(idat).l = 'u-ac043';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284; idat=idat+1;
        % fileUM{idat} = '/u-ab754/u-ab754_Nd_times_LYR_.pp.nc'; labs_UM(idat).l = 'u-ab754';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284; idat=idat+1;

        
end
 

for idat=1:length(fileUM)
    fileUM_rho{idat} = remove_character(fileUM{idat},'VAR_NAME','rho');
    if iadd_umid_label==1
        labs_UM(idat).l = [um_str{idat} ' ' labs_UM(idat).l];
    end
end

vars_out.dirUM = dirUM;
vars_out.fileUM = fileUM;
vars_out.fileUM_rho = fileUM_rho;
vars_out.labs_UM = labs_UM;
vars_out.line_patternUM = line_patternUM;
vars_out.line_colourUM = line_colourUM;
vars_out.marker_styleUM = marker_styleUM;
vars_out.pole_lat = pole_lat;
vars_out.pole_lon = pole_lon;
vars_out.i_aerosol_processing_multi = i_aerosol_processing_multi;
vars_out.append_str_timser = append_str_timser;