%savedir_var='/home/disk/eos1/d.grosvenor/modis_work/saved_data_L3/';
savedir_var='/home/disk/eos8/d.grosvenor/saved_data_L3/';

comp='UWchallenger';

if exist('ioverride_saveload') & ioverride_saveload==1
    clear ioverride_saveload  %use the values set outside and reset the flag

else %otherwise use these values
    saveload='save';
    saveload='load';
end

varsaveload = 'Bootstrap_error_for_TERRA_2008_Nd_with_NP.GTE.50.AND.maxSZA.LT.60.AND.sensZA.LTE.40(days_1-364)'; %actually is 356&357 that are missing

%need to re-do - AQUA data wasn't right
%varsaveload ='Bootstrap_error_for_TERRA+AQUA_2008_Nd_with_NP.GTE.50_AND_CF.GTE.0.8_AND_maxSZA.LT.60_AND_sensZA.LTE.40';
varsaveload = 'Bootstrap_error_for_TERRA+AQUA_2008_Nd_with_NP.GTE.50.AND.sensZA.LTE.40';
varsaveload = 'Bootstrap_error_for_TERRA+AQUA_2008_Nd_with_NP.GTE.50_AND_CF.GTE.0.8_AND_sensZA.LTE.30';

%max_sensZA<=10
varsaveload = 'Bootstrap_error_for_TERRA+AQUA_2008_Nd_with_NP.GTE.50_AND_CF.GTE.0.8_AND_MAXsensZA.LTE.10_20111026T130825.mat';
%max_sensZA<=20
varsaveload = 'Bootstrap_error_for_TERRA+AQUA_2008_Nd_with_NP.GTE.50_AND_CF.GTE.0.8_AND_MAXsensZA.LTE.20_20111026T141616.mat';


varsaveload = ['Bootstrap_error_for_TERRA+AQUA_2008_Nd_with_' thresh_str];

switch saveload
    case 'save'
         savename=[savedir_var varsaveload '_' datestr(now,30) '.mat'];
         savename=remove_character(savename,' ','_');
         save(savename,'boot_out');
    case 'load'
%         switch varsaveload
%             case 'Bootstrap_error_for_TERRA_2008_Nd_with_NP.GTE.50.AND.maxSZA.LT.60.AND.sensZA.LTE.40(days_1-364)';                
%                 filename='Bootstrap_error_for_TERRA_2008_Nd_with_NP.GTE.50.AND.maxSZA.LT.60.AND.sensZA.LTE.40(days_1-364)_20111024T173148';
%                 
%             case 'Bootstrap_error_for_TERRA+AQUA_2008_Nd_with_NP.GTE.50.AND.sensZA.LTE.40' 
%                 filename='Bootstrap_error_for_TERRA+AQUA_2008_Nd_with_NP.GTE.50.AND.sensZA.LTE.40_20111025T103921';
%                 
%             case 'Bootstrap_error_for_TERRA+AQUA_2008_Nd_with_NP.GTE.50_AND_CF.GTE.0.8_AND_maxSZA.LT.60_AND_sensZA.LTE.40'
%                 filename='Bootstrap_error_for_TERRA+AQUA_2008_Nd_with_NP.GTE.50_AND_CF.GTE.0.8_AND_maxSZA.LT.60_AND_sensZA.LTE.40_20111025T133913';
%             case 'Bootstrap_error_for_TERRA+AQUA_2008_Nd_with_NP.GTE.50_AND_CF.GTE.0.8_AND_sensZA.LTE.30'
%                 filename='Bootstrap_error_for_TERRA+AQUA_2008_Nd_with_NP.GTE.50_AND_CF.GTE.0.8_AND_MAXsensZA.LTE.10_20111026T130825.mat';
%         end

        filename = varsaveload;
        
        loadname=[savedir_var filename '.mat'];
        load(loadname);
end

disp('****   Save/Load done   ****');


        
        
     