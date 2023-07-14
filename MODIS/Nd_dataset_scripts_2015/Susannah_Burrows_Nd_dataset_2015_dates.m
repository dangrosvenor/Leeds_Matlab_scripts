%load_file='/home/disk/eos1/d.grosvenor/saved_misc_mat_files/MODIS_multi_year_multi_screening_reff_20150211T232329.mat';
%var_str='reff'; aqua_terra_str='AQUA_AND_TERRA'; year_range=[2000:2014];

%load_file='/home/disk/eos1/d.grosvenor/saved_misc_mat_files/MODIS_multi_year_multi_screening_Nd_20150211T232307.mat';
%var_str='Nd'; aqua_terra_str='AQUA_AND_TERRA'; year_range=[2000:2014];

%load_file='/home/disk/eos1/d.grosvenor/saved_misc_mat_files/MODIS_multi_year_multi_screening_reff_20150211T124750.mat';
%var_str='reff'; aqua_terra_str='AQUA'; year_range=[2002:2014];

%load_file='/home/disk/eos1/d.grosvenor/saved_misc_mat_files/MODIS_multi_year_multi_screening_Nd_20150211T124821.mat';
%var_str='Nd'; aqua_terra_str='AQUA'; year_range=[2002:2014];


%load_file='/home/disk/eos1/d.grosvenor/saved_misc_mat_files/MODIS_multi_year_multi_screening_reff_20150212T014533.mat';
%var_str='reff'; aqua_terra_str='TERRA'; year_range=[2000:2014];

%load_file='/home/disk/eos1/d.grosvenor/saved_misc_mat_files/MODIS_multi_year_multi_screening_Nd_20150212T014544.mat';
%var_str='Nd'; aqua_terra_str='TERRA'; year_range=[2000:2014];

%load_file='/home/disk/eos1/d.grosvenor/saved_misc_mat_files/MODIS_multi_year_multi_screening_reff_20150216T071222.mat';
%var_str='reff'; aqua_terra_str='TERRA'; year_range=[2000:2014];

load_file='/home/disk/eos1/d.grosvenor/saved_misc_mat_files/MODIS_multi_year_multi_screening_Nd_20150216T071248.mat';
var_str='Nd'; aqua_terra_str='TERRA'; year_range=[2000:2014]; save_case='Florent';

%all CF, CTT > 268 (for the CTT screenings 
load_file='/home/disk/eos1/d.grosvenor/saved_misc_mat_files/MODIS_multi_year_multi_screening_Nd_20150402T131117.mat';
var_str='Nd'; aqua_terra_str='TERRA'; year_range=[2000:2014]; save_case='UKESM_Jane';
screenings = {'no_screening','CTT_gt_268','SZA_lt_70','SZA_lt_70_Nd_AND_CTT_gt_268'};

%CF>80, CTT>273.15
load_file='/home/disk/eos1/d.grosvenor/saved_misc_mat_files/MODIS_multi_year_multi_screening_Nd_20150616T052531.mat'; 
var_str='Nd'; aqua_terra_str='TERRA'; year_range=[2000:2014]; save_case='UKESM_Jane';
%Specify the screenings that relate to ones in the file
screenings = {'no_screening','CTT_gt_273_CF_gt_80','SZA_lt_70_CF_gt_80','SZA_lt_70_CF_gt_80_AND_CTT_gt_273'};

% all CF, CTT>273.15
load_file='/home/disk/eos1/d.grosvenor/saved_misc_mat_files/MODIS_multi_year_multi_screening_Nd_20150618T062603.mat'; 
var_str='Nd'; aqua_terra_str='TERRA'; year_range=[2000:2014]; save_case='UKESM_Jane';
%Specify the screenings that relate to ones in the file
screenings = {'no_screening','CTT_gt_273_allCF','SZA_lt_70_allCF','SZA_lt_70_allCF_AND_CTT_gt_273'};

% Jan 2016 - production of Nd dataset for Florent and Jim since want to use
% Aqua instead of Terra.
load_file='/home/disk/eos1/d.grosvenor/saved_misc_mat_files/MODIS_multi_year_multi_screening_Nd_20160108T033901.mat'; %Test file with just 2002
%var_str='Nd'; aqua_terra_str='AQUA'; year_range=[2000:2014]; save_case='Florent';
var_str='Nd'; aqua_terra_str='AQUA'; year_range=[2002:2014]; save_case='Florent';

fprintf(1,'\nSaving...');


load(load_file);


savefile=['/home/disk/eos1/d.grosvenor/saved_misc_mat_files/MODIS_' aqua_terra_str '_' var_str '_QA_from_L3_daily_' save_case '_2000_2014_' datestr(now,30)];
savefile2=[savefile '.mat'];
        

i=0;
for iyear=year_range
    i=i+1;
    year_str = ['y' num2str(iyear)];   
    
    for j=1:length(eval(['Psave_' year_str]))
        Psave{j}(:,:,i) = eval(['Psave_' year_str '{j}(:,:)']);
        Nsave{j}(:,:,i) = eval(['Nsave_' year_str '{j}(:,:)']);        
        Nd_std_dev_save{j}(:,:,i) = eval(['Nd_std_dev_save_' year_str '{j}(:,:)']);            
    end
    
    
end



YEARS = year_range;

save(savefile2,'MLAT','MLON','YEARS');

%function save_vars_for_Florent(savefile,Psave,Nsave,ind,varname)
  %ind is the required index for the multi-screening cell
  
switch save_case
    case 'Florent'
        
        
        %data is ordered by month first and then the different screenings.
        %So here there were 3 months and 9 screenings - 27 in total.
        %Screenings 2-5 should be a repeat of 6-9 (but thresh values different). So start at screening 6,
        %which equates to no.16 onwards (6-1)*3 + 1 = 16
        
        %This writes out the variables for each month separately. Perhaps
        %better to put all the months for each year together?? Actually
        %maybe not since this makes doing e.g. JJA easier.
  
save_vars_for_Florent_func(savefile2,Psave,Nsave,16,['Aug_no_screening_' var_str]);    
save_vars_for_Florent_func(savefile2,Psave,Nsave,17,['Sep_no_screening_' var_str]);    
save_vars_for_Florent_func(savefile2,Psave,Nsave,18,['Oct_no_screening_' var_str]);    

save_vars_for_Florent_func(savefile2,Psave,Nsave,19,['Aug_CTT_gt_268_' var_str]);    
save_vars_for_Florent_func(savefile2,Psave,Nsave,20,['Sep_CTT_gt_268_' var_str]);    
save_vars_for_Florent_func(savefile2,Psave,Nsave,21,['Oct_CTT_gt_268_' var_str]);    

save_vars_for_Florent_func(savefile2,Psave,Nsave,22,['Aug_SZA_lt_70_' var_str]);    
save_vars_for_Florent_func(savefile2,Psave,Nsave,23,['Sep_SZA_lt_70_' var_str]);    
save_vars_for_Florent_func(savefile2,Psave,Nsave,24,['Oct_SZA_lt_70_' var_str]);    

save_vars_for_Florent_func(savefile2,Psave,Nsave,25,['Aug_SZA_lt_70_Nd_AND_CTT_gt_268_' var_str]);    
save_vars_for_Florent_func(savefile2,Psave,Nsave,26,['Sep_SZA_lt_70_Nd_AND_CTT_gt_268_' var_str]);    
save_vars_for_Florent_func(savefile2,Psave,Nsave,27,['Oct_SZA_lt_70_Nd_AND_CTT_gt_268_' var_str]);  

    case 'UKESM_Jane'        
        
        months={'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};

        nscreen=9;
        %indices of the screenings relative to the total number (nscreen):-
        iscreen=6:9;
        for is=1:length(screenings)
            for im=1:length(months)
                month=months{im};
                screen=screenings{is};
                is2=iscreen(is);
                is3 = (is2-1)*length(months)+im
                save_vars_for_Florent_func2(savefile2,Psave,Nsave,Nd_std_dev_save,is3,[month '_' screen '_' var_str]); 


            end
        end
        
        
        
end


mat2nc_Dan(savefile2,[savefile '.nc']);




fprintf(1,' Done\n');
