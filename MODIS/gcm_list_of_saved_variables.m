clear gcm_var; istring=1;

if  load_gcm_process_vars_flag==0 & save_gcm_process_vars_flag==0
    textinput = input('Don''t forget to set the time_dim value!! Not loading/saving any gcm_processed variables. Press enter to continue');
end


%LTS related
%gcm_var{istring}='gcm_theta';  field_string{istring}='.timeseries3'; time_dim{istring}=1; istring=istring+1;%
gcm_var{istring}='gcm_theta700';  field_string{istring}='.timeseries3'; time_dim{istring}=1; istring=istring+1;%
gcm_var{istring}='gcm_theta0';  field_string{istring}='.timeseries3'; time_dim{istring}=1; istring=istring+1;%
gcm_var{istring}='gcm_LTS';  field_string{istring}='.timeseries3'; time_dim{istring}=1; istring=istring+1;%
gcm_var{istring}='gcm_qv700';  field_string{istring}='.timeseries3'; time_dim{istring}=1; istring=istring+1;%
%gcm_var{istring}='gcm_theta1000';  field_string{istring}='.timeseries3'; time_dim{istring}=1; istring=istring+1;%
%gcm_var{istring}='gcm_LTS1000';  field_string{istring}='.timeseries3'; time_dim{istring}=1; istring=istring+1;%


%  gcm_var{istring}='rwp';  field_string{istring}='.timeseries3'; time_dim{istring}=1; istring=istring+1;%

gcm_var{istring}='gcm_time_matlab'; field_string{istring}='.timeseries3'; time_dim{istring}=1; istring=istring+1;%
gcm_var{istring}='gcm_time_days'; field_string{istring}='.timeseries3'; time_dim{istring}=1; istring=istring+1;%
gcm_var{istring}='gcm_time_UTC'; field_string{istring}='.timeseries3'; time_dim{istring}=1; istring=istring+1;%
gcm_var{istring}='gcm_month'; field_string{istring}='.timeseries3'; time_dim{istring}=1; istring=istring+1;%
gcm_var{istring}='gcm_decimal_days'; field_string{istring}='.timeseries3'; time_dim{istring}=1; istring=istring+1;%
gcm_var{istring}='daynum_timeseries3'; field_string{istring}='.timeseries3'; time_dim{istring}=1; istring=istring+1;%

gcm_var{istring}='gcm_Plat2D'; field_string{istring}='.timeseries3'; time_dim{istring}=0; istring=istring+1;%
gcm_var{istring}='gcm_Plon2D'; field_string{istring}='.timeseries3'; time_dim{istring}=0; istring=istring+1;%
gcm_var{istring}='gcm_Plat2D_edges'; field_string{istring}='.timeseries3'; time_dim{istring}=0; istring=istring+1;%
gcm_var{istring}='gcm_Plon2D_edges'; field_string{istring}='.timeseries3'; time_dim{istring}=0; istring=istring+1;%

gcm_var{istring}='gcm_lwp';  field_string{istring}='.timeseries3'; time_dim{istring}=1; istring=istring+1;%

%gcm_var{istring}='gcm_strat_Nd';  field_string{istring}='.timeseries3'; time_dim{istring}=1; istring=istring+1;%
%gcm_var{istring}='gcm_REFFL';  field_string{istring}='.timeseries3'; time_dim{istring}=1; istring=istring+1;%


        gcm_var{istring}='iwp_isccp_low';  field_string{istring}='.timeseries3'; time_dim{istring}=1; istring=istring+1;%
        gcm_var{istring}='lwp_isccp_low';  field_string{istring}='.timeseries3'; time_dim{istring}=1; istring=istring+1;%

        gcm_var{istring}='iwp_isccp_mid';  field_string{istring}='.timeseries3'; time_dim{istring}=1; istring=istring+1;%
        gcm_var{istring}='lwp_isccp_mid';  field_string{istring}='.timeseries3'; time_dim{istring}=1; istring=istring+1;%

        gcm_var{istring}='iwp_isccp_high';  field_string{istring}='.timeseries3'; time_dim{istring}=1; istring=istring+1;%
        gcm_var{istring}='lwp_isccp_high';  field_string{istring}='.timeseries3'; time_dim{istring}=1; istring=istring+1;%
         
        
    gcm_var{istring}='gcm_landmask';  field_string{istring}='.timeseries3'; time_dim{istring}=0; istring=istring+1;%
    
            %RWP related (needs 3D rain) - only have for CAM5 and CAMCLUBBv2 at the moment
        gcm_var{istring}='rwp';  field_string{istring}='.timeseries3'; time_dim{istring}=1; istring=istring+1;%
        gcm_var{istring}='rwp_isccp_low';  field_string{istring}='.timeseries3'; time_dim{istring}=1; istring=istring+1;%
        gcm_var{istring}='rwp_isccp_mid';  field_string{istring}='.timeseries3'; time_dim{istring}=1; istring=istring+1;%        
        gcm_var{istring}='rwp_isccp_high';  field_string{istring}='.timeseries3'; time_dim{istring}=1; istring=istring+1;%        
   
        


%for load_gcm_process_vars only load the bare min to save memory
if load_gcm_process_vars_flag~=4

    gcm_var{istring}='gcm_landmask';  field_string{istring}='.timeseries3'; time_dim{istring}=0; istring=istring+1;%
    gcm_var{istring}='gcm_zsurf';  field_string{istring}='.timeseries3'; time_dim{istring}=0; istring=istring+1;%

    gcm_var{istring}='gcm_ps';  field_string{istring}='.timeseries3'; time_dim{istring}=1; istring=istring+1;%%- = 1e-6 * gcm_drop; %cm3 - already done the /kg to /m3 conversion
    gcm_var{istring}='gcm_tsurf';  field_string{istring}='.timeseries3'; time_dim{istring}=1; istring=istring+1;%


    gcm_var{istring}='gcm_iwp';  field_string{istring}='.timeseries3'; time_dim{istring}=1; istring=istring+1;%
    %gcm_var{istring}='gcm_lwp_all_clouds';  field_string{istring}='.timeseries3'; time_dim{istring}=1; istring=istring+1;%
    %gcm_var{istring}='gcm_iwp_all_clouds';  field_string{istring}='.timeseries3'; time_dim{istring}=1; istring=istring+1;%
    %what is gcm_lwp_all_clouds??

 
end



%always save these, but only load in certain cases (3D/4D fields)
if save_gcm_process_vars_flag>0 |  load_gcm_process_vars_flag==2 
    gcm_var{istring}='gcm_phalf';  field_string{istring}='.timeseries3'; time_dim{istring}=1; istring=istring+1;%
    gcm_var{istring}='gcm_pfull';  field_string{istring}='.timeseries3'; time_dim{istring}=1; istring=istring+1;%
    gcm_var{istring}='gcm_temp';  field_string{istring}='.timeseries3'; time_dim{istring}=1; istring=istring+1;%
    gcm_var{istring}='gcm_rho';  field_string{istring}='.timeseries3'; time_dim{istring}=1; istring=istring+1;%    
    gcm_var{istring}='gcm_rain3D';  field_string{istring}='.timeseries3'; time_dim{istring}=1; istring=istring+1;%    
    gcm_var{istring}='gcm_drop_read';   field_string{istring}='.timeseries3'; time_dim{istring}=1; istring=istring+1;%             
    gcm_var{istring}='gcm_cf';  field_string{istring}='.timeseries3'; time_dim{istring}=1; istring=istring+1;%
    
 %3D variables -  if load_gcm_process_vars_flag~=4
        gcm_var{istring}='gcm_liq_av';  field_string{istring}='.timeseries3'; time_dim{istring}=1; istring=istring+1;%
        gcm_var{istring}='gcm_drop2';   field_string{istring}='.timeseries3'; time_dim{istring}=1; istring=istring+1;%
        
        
end


     


        

if save_gcm_process_vars_flag>0
   %can move ones that don't want to load, but want to save in the future
   %in here - put the loading flag in a comment after
   %gcm_var{istring}='gcm_cloudsat_CFAD'; e.g.:-
   %field_string{istring}='.timeseries3'; time_dim{istring}=1; istring=istring+1;%if save_gcm_process_vars_flag>0 | cosp_flag4D==1

   %radiation fields
gcm_var{istring}='LW_surf_down';  field_string{istring}='.timeseries3'; time_dim{istring}=1; istring=istring+1;%
gcm_var{istring}='LW_surf_net';  field_string{istring}='.timeseries3'; time_dim{istring}=1; istring=istring+1;%
gcm_var{istring}='LW_TOA_net';  field_string{istring}='.timeseries3'; time_dim{istring}=1; istring=istring+1;%
gcm_var{istring}='SW_surf_down';  field_string{istring}='.timeseries3'; time_dim{istring}=1; istring=istring+1;%
gcm_var{istring}='SW_surf_net';  field_string{istring}='.timeseries3'; time_dim{istring}=1; istring=istring+1;%
gcm_var{istring}='SW_TOA_net';  field_string{istring}='.timeseries3'; time_dim{istring}=1; istring=istring+1;%
gcm_var{istring}='SW_TOA_up';  field_string{istring}='.timeseries3'; time_dim{istring}=1; istring=istring+1;%
gcm_var{istring}='albedo';  field_string{istring}='.timeseries3'; time_dim{istring}=1; istring=istring+1;%

  


                 
  
   

   


       
        

        gcm_var{istring}='gcm_REFFL_max_noCF';  field_string{istring}='.timeseries3'; time_dim{istring}=1; istring=istring+1;%
        gcm_var{istring}='gcm_REFFL_maxliq';  field_string{istring}='.timeseries3'; time_dim{istring}=1; istring=istring+1;%
        gcm_var{istring}='gcm_REFFL_maxlayer';  field_string{istring}='.timeseries3'; time_dim{istring}=1; istring=istring+1;%
        
end

 %these had no conditions for loading


      

if load_gcm_process_vars_flag~=4

end

    
    if ((exist('iread_cosp') & iread_cosp==1) & save_gcm_process_vars_flag>0) | cosp_flag4D==1
        if ~exist('ino_dbz_cfads') | ino_dbz_cfads==0
            gcm_var{istring}='gcm_cloudsat_CFAD';  field_string{istring}='.timeseries3'; time_dim{istring}=1; istring=istring+1;%
            gcm_var{istring}='gcm_calipso_CFAD';  field_string{istring}='.timeseries3'; time_dim{istring}=1; istring=istring+1;%
        end
        

        if ~exist('ino_dbz_cfads') | ino_dbz_cfads==0
            gcm_var{istring}='CFAD_dbz_edges';  field_string{istring}='.timeseries3'; time_dim{istring}=0; istring=istring+1;%
            gcm_var{istring}='cf_CFAD_dbz';  field_string{istring}='.timeseries3'; time_dim{istring}=1; istring=istring+1;%
            gcm_var{istring}='mean_CFAD_dbz';  field_string{istring}='.timeseries3'; time_dim{istring}=1; istring=istring+1;%
            gcm_var{istring}='dbz_thresh';  field_string{istring}='.timeseries3'; time_dim{istring}=0; istring=istring+1;%
        
        
        gcm_var{istring}='cfad_alts_edges';  field_string{istring}='.timeseries3'; time_dim{istring}=0; istring=istring+1;%
        gcm_var{istring}='CFAD_sr_edges';  field_string{istring}='.timeseries3'; time_dim{istring}=0; istring=istring+1;%
        
       
        gcm_var{istring}='cf_CFAD_sr';  field_string{istring}='.timeseries3'; time_dim{istring}=1; istring=istring+1;%
        gcm_var{istring}='mean_CFAD_sr';  field_string{istring}='.timeseries3'; time_dim{istring}=1; istring=istring+1;%
        gcm_var{istring}='sr_thresh';  field_string{istring}='.timeseries3'; time_dim{istring}=0; istring=istring+1;%
      
        end

    end


%if save_gcm_process_vars_flag==1 | load_gcm_process_vars_flag==1 | load_gcm_process_vars_flag==2  %3D fields that are slow to calculate - perhaps won't be produced    
if save_gcm_process_vars_flag==1 | load_gcm_process_vars_flag==1  %3D fields that are slow to calculate - perhaps won't be produced        
        gcm_var{istring}='h_half';  field_string{istring}='.timeseries3'; time_dim{istring}=1; istring=istring+1;%
        gcm_var{istring}='h_full';  field_string{istring}='.timeseries3'; time_dim{istring}=1; istring=istring+1;%
end  

if cosp_flag==1      
    gcm_var{istring}='cllcalipso';  field_string{istring}='.timeseries3'; time_dim{istring}=1; istring=istring+1;%
    gcm_var{istring}='clmcalipso';  field_string{istring}='.timeseries3'; time_dim{istring}=1; istring=istring+1;%
    gcm_var{istring}='clhcalipso';  field_string{istring}='.timeseries3'; time_dim{istring}=1; istring=istring+1;%


    gcm_var{istring}='liqCF_modis';  field_string{istring}='.timeseries3'; time_dim{istring}=1; istring=istring+1;%
    gcm_var{istring}='liqRe_modis';  field_string{istring}='.timeseries3'; time_dim{istring}=1; istring=istring+1;%
    gcm_var{istring}='liqTau_modis';  field_string{istring}='.timeseries3'; time_dim{istring}=1; istring=istring+1;%    
end



%fields from gcm_process
if save_gcm_process_vars_flag>0 | load_gcm_process_vars_flag>0

    
    if load_gcm_process_vars_flag~=4
        
      
    end
    
        gcm_var{istring}='gcm_precL';  field_string{istring}='.timeseries3'; time_dim{istring}=1; istring=istring+1;%
        gcm_var{istring}='gcm_precT';  field_string{istring}='.timeseries3'; time_dim{istring}=1; istring=istring+1;%
        


        gcm_var{istring}='gcm_Nd_max_noCF';  field_string{istring}='.timeseries3'; time_dim{istring}=1; istring=istring+1;%
        gcm_var{istring}='gcm_Nd_maxliq';  field_string{istring}='.timeseries3'; time_dim{istring}=1; istring=istring+1;%                 %- Nd at max LWC over all levels (after screening)
        gcm_var{istring}='gcm_Nd_meanlayer';  field_string{istring}='.timeseries3'; time_dim{istring}=1; istring=istring+1;%              %  - mean Nd over the the identified continuous liquid cloud layer
        gcm_var{istring}='gcm_Nd_max_lwc_cont';  field_string{istring}='.timeseries3'; time_dim{istring}=1; istring=istring+1;%           %- Nd at the location of the max LWC within the continuous cloud layer
        gcm_var{istring}='gcm_LWC_max';  field_string{istring}='.timeseries3'; time_dim{istring}=1; istring=istring+1;%
        gcm_var{istring}='gcm_CF_maxliq';  field_string{istring}='.timeseries3'; time_dim{istring}=1; istring=istring+1;%
        gcm_var{istring}='gcm_CF_max_screened';  field_string{istring}='.timeseries3'; time_dim{istring}=1; istring=istring+1;%

        gcm_var{istring}='gcm_CTT_layer';  field_string{istring}='.timeseries3'; time_dim{istring}=1; istring=istring+1;%

        gcm_var{istring}='gcm_CTP';  field_string{istring}='.timeseries3'; time_dim{istring}=1; istring=istring+1;%
        gcm_var{istring}='gcm_CBP';  field_string{istring}='.timeseries3'; time_dim{istring}=1; istring=istring+1;%






    gcm_var{istring}='gcm_Nd_maxlayer';  field_string{istring}='.timeseries3'; time_dim{istring}=1; istring=istring+1;%               %     - max Nd in the continuous layer
    gcm_var{istring}='gcm_lwp_minthreshCF';  field_string{istring}='.timeseries3'; time_dim{istring}=1; istring=istring+1;%
    gcm_var{istring}='gcm_Nd_max_screen';  field_string{istring}='.timeseries3'; time_dim{istring}=1; istring=istring+1;%             %- max Nd in all layers (after screening for whatever screening was requested from the switches
    
%    gcm_var{istring}='gcm_Nd_max_screen_isccp_low';  field_string{istring}='.timeseries3'; time_dim{istring}=1; istring=istring+1;%             %- max Nd in all layers (after screening for whatever screening was requested from the switches
%    gcm_var{istring}='gcm_Nd_max_screen_isccp_mid';  field_string{istring}='.timeseries3'; time_dim{istring}=1; istring=istring+1;%             %- max Nd in all layers (after screening for whatever screening was requested from the switches
%    gcm_var{istring}='gcm_Nd_max_screen_isccp_high';  field_string{istring}='.timeseries3'; time_dim{istring}=1; istring=istring+1;%             %- max Nd in all layers (after screening for whatever screening was requested from the switches

%    gcm_var{istring}='gcm_Nd_max_noCF_isccp_low';  field_string{istring}='.timeseries3'; time_dim{istring}=1; istring=istring+1;%
%    gcm_var{istring}='gcm_Nd_max_noCF_isccp_mid';  field_string{istring}='.timeseries3'; time_dim{istring}=1; istring=istring+1;%
%    gcm_var{istring}='gcm_Nd_max_noCF_isccp_high';  field_string{istring}='.timeseries3'; time_dim{istring}=1; istring=istring+1;%
    

    %isccp style cloud fractions broken into low, mid and high cloud, along
    %with iwps calculated over these ranges
    gcm_var{istring}='cf_isccp_low';  field_string{istring}='.timeseries3'; time_dim{istring}=1; istring=istring+1;%
    gcm_var{istring}='cf_isccp_mid';  field_string{istring}='.timeseries3'; time_dim{istring}=1; istring=istring+1;%
    gcm_var{istring}='cf_isccp_high';  field_string{istring}='.timeseries3'; time_dim{istring}=1; istring=istring+1;%

    %before and after mphys LWPs
    if ilwcAPBP==1
        gcm_var{istring}='lwpBP_isccp_low';  field_string{istring}='.timeseries3'; time_dim{istring}=1; istring=istring+1;%
        gcm_var{istring}='lwpAP_isccp_low';  field_string{istring}='.timeseries3'; time_dim{istring}=1; istring=istring+1;%
        gcm_var{istring}='lwpBP_isccp_mid';  field_string{istring}='.timeseries3'; time_dim{istring}=1; istring=istring+1;%
        gcm_var{istring}='lwpAP_isccp_mid';  field_string{istring}='.timeseries3'; time_dim{istring}=1; istring=istring+1;%
        gcm_var{istring}='lwpBP_isccp_high';  field_string{istring}='.timeseries3'; time_dim{istring}=1; istring=istring+1;%
        gcm_var{istring}='lwpAP_isccp_high';  field_string{istring}='.timeseries3'; time_dim{istring}=1; istring=istring+1;%                
        gcm_var{istring}='lwpSEDTEN_isccp_low';  field_string{istring}='.timeseries3'; time_dim{istring}=1; istring=istring+1;%
        gcm_var{istring}='lwpSEDTEN_isccp_mid';  field_string{istring}='.timeseries3'; time_dim{istring}=1; istring=istring+1;%
        gcm_var{istring}='lwpSEDTEN_isccp_high';  field_string{istring}='.timeseries3'; time_dim{istring}=1; istring=istring+1;%
        gcm_var{istring}='lwpEVAPTEN_isccp_low';  field_string{istring}='.timeseries3'; time_dim{istring}=1; istring=istring+1;%
        gcm_var{istring}='lwpEVAPTEN_isccp_mid';  field_string{istring}='.timeseries3'; time_dim{istring}=1; istring=istring+1;%
        gcm_var{istring}='lwpEVAPTEN_isccp_high';  field_string{istring}='.timeseries3'; time_dim{istring}=1; istring=istring+1;%
    end
    
    if save_gcm_process_vars_flag==1 | load_gcm_process_vars_flag==1  %sometimes don't have the fields that require height
        %can an estimate of LWP be integrated over pressure instead? Yes,
        %see gcm_process! (dp/dz = -rho*g, so can calc dz)

       

        gcm_var{istring}='gcm_CTH';  field_string{istring}='.timeseries3'; time_dim{istring}=1; istring=istring+1;%
        gcm_var{istring}='gcm_CBH';  field_string{istring}='.timeseries3'; time_dim{istring}=1; istring=istring+1;%


    end
    
       
            

        

end



  


%clear gcm_var; istring=1;

