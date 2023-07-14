
% Runs ACSIS_Robson_paper_load_data for each variable and then runs 
% ACSIS_Robson_paper_plot_timeseries_RUN_multi_generic

%MIP = 'AMIP';
MIPs = {'DAMIP','CMIP','CMIP with obs','UKESM1-AerChemMIP_control'};
MIPs = {'CMIP with obs','UKESM1-AerChemMIP_control'};
MIPs = {'CMIP with obs'};
%MIPs = {'UKESM1-AerChemMIP_control'};
%MIPs = {'HADGEM-LL'};
%MIPs = {'DAMIP'};
%MIPs = {'CMIP'}; %(UKESM full ensemble)
%MIPs = {'UKESM1-AerChemMIP_control'}; %Also for pi-aer and aerosol-only proxy
%MIPs = {'HADGEM-LL'};

%vars_ukesm = {'rsutcs','scldncl','rsut','lwp','lwpic','clt','ts',};
vars_ukesm = {'rsutcs','scldncl','rsut','lwpic','clt','ts',};
vars_ukesm = {'rsut','rsutcs'};
vars_ukesm = {'rsutcs'};
%vars_ukesm = {'scldncl'};
vars_ukesm = {'rsutcs'};
vars_ukesm = {'clt'};
vars_ukesm = {'lwpic'};
vars_ukesm = {'lwp'};
%vars_ukesm = {'ts'};
%vars_ukesm = {'rsutcs','lwpic','clt','ts'};
%vars_ukesm = {'rsut'};
%vars_ukesm = {'Nd_clw_weighted_ESGF'};
%vars_ukesm = {'od550aer'};
%vars_ukesm = {'od550tot'};
vars_ukesm = {'rsut','Nd_clw_weighted_ESGF','od550tot','clt','lwp','ts',}; %ones with obs for paper

for iMIP=1:length(MIPs)
    MIP_temp = MIPs{iMIP}; 
    MIP = MIP_temp;
    switch MIP_temp
        %case {'CMIP','HADGEM-LL'}
            %iadd_DAMIP = 0;
            %MIP = MIP_temp;
            %obs_str_DRIVER = '1850-1970 1971-2014 plot';
            %obs_str = 'default'; %Should plot the 3 periods
        case 'DAMIP'
            iadd_DAMIP = 1;
            MIP = 'HADGEM-LL';
            %obs_str_DRIVER = '1850-1970 1971-2014 no plot';
       case 'CMIP with obs'            
            MIP = 'CMIP';
            %obs_str_DRIVER = '1850-1970 1971-2014 no plot';     
    end

    for iv=1:length(vars_ukesm)                
        var_ukesm = vars_ukesm{iv};
        
        
% --------------- LOAD the DATA ------------------------------        
        %Run loading script
        ioverride_vals=1;
        ACSIS_Robson_paper_load_data
        
        switch MIP_temp
            case 'UKESM1-AerChemMIP_control'
                iload_DAMIP = 0; %only actually does this if iadd_DAMIP==1 too
                iadd_DAMIP = 0;
                
                iload_AerChemMIP=1; %only actually does this if iadd_AerChemMIP=1 too
                iadd_AerChemMIP=1;
                
                i_CMIP6_multimodel=0;
                iload_multi_model=0;
                
                iadd_HADGEM=0; %
                iload_HADGEM=0;
                
                iadd_amip=0;
                iload_amip=0;
                
                iadd_nudged = 0;
                
            case 'DAMIP'
                iload_DAMIP = 1; %only actually does this if iadd_DAMIP==1 too
                iadd_DAMIP = 1;
                
                iload_AerChemMIP=0; %only actually does this if iadd_AerChemMIP=1 too
                iadd_AerChemMIP=0;
                
                i_CMIP6_multimodel=0;
                iload_multi_model=0;
                
                iadd_HADGEM=0; %
                iload_HADGEM=0;
                
                iadd_amip=0;
                iload_amip=0;
                
                iadd_nudged = 0;
                
                
            case {'CMIP','CMIP with obs'}
                                
                iload_DAMIP = 0; %only actually does this if iadd_DAMIP==1 too
                iadd_DAMIP = 0;
                
                iload_AerChemMIP=0; %only actually does this if iadd_AerChemMIP=1 too
                iadd_AerChemMIP=0;
                
                i_CMIP6_multimodel=0;
                iload_multi_model=0;
                
                iadd_HADGEM=1; %
                iload_HADGEM=1;
                
                iadd_amip=0;
                iload_amip=0;
                
                iadd_nudged = 0;
                
                switch MIP_temp
                    case 'CMIP'
                        obs_str='';
                    case 'CMIP with obs'
                        switch var_ukesm
                            case {'rsut'}
                                iadd_amip=1;
                                iload_amip=1;                                
                        end
   
                        
                end
                
            otherwise 
                                
                iload_DAMIP = 0; %only actually does this if iadd_DAMIP==1 too
                iadd_DAMIP = 0;
                
                iload_AerChemMIP=0; %only actually does this if iadd_AerChemMIP=1 too
                iadd_AerChemMIP=0;
                
                i_CMIP6_multimodel=0;
                iload_multi_model=0;
                
                iadd_HADGEM=0; %
                iload_HADGEM=0;
                
                iadd_amip=0;
                iload_amip=0;
                
                iadd_nudged = 0;
                
                
        end
                
        
% --------------- Run the timeseries scripts ------------------------------

        %obs_str='none';
        %obs_str = obs_str_DRIVER;
        ioverride_vals=1;
        ACSIS_Robson_paper_plot_timeseries_RUN_multi_generic
        
%         switch MIP
%             case 'UKESM1-AerChemMIP_control'
%                 iload_DAMIP = 0; %only actually does this if iadd_DAMIP==1 too
%                 iadd_DAMIP = 0;
%                 
%                 iload_AerChemMIP=1; %only actually does this if iadd_AerChemMIP=1 too
%                 iadd_AerChemMIP=1;
%                 
%                 i_CMIP6_multimodel=0;
%                 iload_multi_model=0;
%                 
%                 iadd_HADGEM=0; %
%                 iload_HADGEM=0;
%                 
%                 ioverride_vals=1;
%                 ACSIS_Robson_paper_plot_timeseries_RUN_multi_generic
%         end
        
        
    end
    
end



% 
% var_ukesm = 'Nd_cf_weighted_UKESM';
% var_ukesm = 'Nd_clw_weighted_ESGF';
% var_ukesm = 'Nd_clw_weighted_ESGF_no_dz';
% %var_ukesm = 'calipso_low_cloud_amount';
% %var_ukesm = 'calipso_total_cloud_amount';
% %var_ukesm = 'SW_up_TOA';
% %var_ukesm = 'DEEPC_fluxes'; %SW TOA fluxes from the DeepC dataset.
% %var_ukesm = 'SO2_low_anthropogenic_emissions';
% %var_ukesm = 'clt';
% var_ukesm = 'rsut'; %Upwelling TOA SW
% %var_ukesm = 'rsutcs'; %Clear-sky upwelling TOA SW
% %var_ukesm = 'rsds'; %Downwelling surface SW
% %var_ukesm = 'ts'; %surface temperature
% %var_ukesm = 'dust_od550';
% %var_ukesm = 'scldncl'; %cloud top stratiform Nd averaged from daily to monthly by me (ACSIS_Robson_paper_calc_monthly_Nd_ESGF.m)
% %N.B. - didn't do liq cloud frac weighting for time average here - not
% %really needed if consider that cumulus should weight equally to
% %stratocumlus
% %var_ukesm='clwvi'; clwvi includes IWP too.
% %var_ukesm='lwp'; %
% %var_ukesm='lwpic'; % N.B. - have divided by the total cloud fraction for this variable (in ACSIS_dat_trends_load_ensemble_esgf.m)
% % to approximate the in-cloud LWP.
% %var_ukesm='prw'; %column water vapour (water vapour path)
% %var_ukesm='od550aer'; %Optical depth at 550nm
% %var_ukesm='od550aerso'; %Optical depth at 550nm
