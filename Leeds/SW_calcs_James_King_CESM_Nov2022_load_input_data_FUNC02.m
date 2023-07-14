function [output] = SW_calcs_James_King_CESM_Nov2022_load_input_data_FUNC02(UM_run,file_prefix,file_suffix,optionals)

%Called from SW_calcs_James_King_CESM_Nov2022_load_input_data.m

%Convert all of the variable names in the input structure to actual names
%for ease of use
name_struc='optionals'; %The name of the structure
names = eval(['fieldnames(' name_struc ');']);
for i=1:length(names)
    eval_str = [names{i} ' = ' name_struc '.' names{i} ';'];
    eval(eval_str);
end

if ~exist('iget_time_only')
    iget_time_only = 0;
end
if ~exist('file_suffix')
    file_suffix = ''; %e.g., '_yrs3-5' for the limited data files
end


% Called from SW_calcs_James_King_CESM_Nov2022_load_input_data.m

data_dir = '/home/disk/eos15/d.grosvenor/UM/James_Weber/James_King_CESM/';
%file_prefix = 'CAM_chem_SSP3-7.0_2050_maxforest_complete_adjust_SP_mode.';


%% 1) Load ensemble monthly mean data
%if iNd_only==0
    %CAM_chem_SSP3-7.0_2050_maxforest_complete_adjust_SP_mode.002_FREQL_yrs3-5.nc
    %CAM_chem_SSP3-7.0_2095_maxforest_complete_adjust_SP_mode.001_LWC_yrs3-5.nc
    
    
    % -- SW TOA up (to compare against the calculation) --    
    %stash = ''; var='UM_m01s01i217_vn1105';
    %SW_calcs_James_Weber_Jan2022_load_input_data_FUNC
    
    var='FSUTOA'; stash = ['002_' var file_suffix];    
    SW_calcs_James_King_CESM_Nov2022_load_input_data_FUNC      
    
    output.nT = size(dat,1);
    if ~exist('nT')
        nT = output.nT;
    end
    
    output.SW_up_TOA = dat(1:nT,:,:); %[25x192x288 double]
    
if iget_time_only == 0
    
    %SW down TOA - calculate from the net and upwelling (net = down minus
    %up, so down = net plus up
    var='FSNTOA'; stash = ['002_' var file_suffix];    
    SW_calcs_James_King_CESM_Nov2022_load_input_data_FUNC  
    output.SW_down_TOA = dat(1:nT,:,:) + output.SW_up_TOA; %[25x192x288 double]
       
    
    
    
    % --- SW TOA upwelling CLEAR-SKY - for use directly for overall SWup calculations : (1-fc) part--- 
    % Only have the net clear-sky TOA. Downwelling shoudl be the same for
    % clear-sky as for all-sky, so use this to calc upwelling.
    % up = down minus net
    var='FSNTOAC'; stash = ['002_' var file_suffix];    
    SW_calcs_James_King_CESM_Nov2022_load_input_data_FUNC  
    output.SW_up_TOA_clear_sky= output.SW_down_TOA - dat(1:nT,:,:) ; %[25x192x288 double]  
            
    % --- SW surface downwelling CLEAR-SKY - for calculation of surface albedo (not used for transmissivity in this setup) ---    
    var='FSDSC'; stash = ['002_' var file_suffix];    
    SW_calcs_James_King_CESM_Nov2022_load_input_data_FUNC  
    output.SW_down_surf_clear_sky = dat(1:nT,:,:);    
    
%     % --- SW surface upwelling CLEAR-SKY - for calculation of surface albedo (not used for transmissivity in this setup) ---
%     % Calculate from downwelling and net (presumably) downwelling) net =
%     % down minus up. up = down minus net
%     var='FSNSC'; stash = ['002_' var file_suffix];    
%     SW_calcs_James_King_CESM_Nov2022_load_input_data_FUNC  
%     output.SW_up_surf_clear_sky = output.SW_down_surf_clear_sky - dat(1:nT,:,:) ; %[25x192x288 double] 
        
    output.SW_up_surf_clear_sky = output.SW_down_surf_clear_sky; %dummy values for now
    
    
    
    
    % --- total CF ---
    var='CLDTOT'; stash = ['002_' var file_suffix];    
    SW_calcs_James_King_CESM_Nov2022_load_input_data_FUNC  
    output.totCF = dat(1:nT,:,:);
    
    
   
    
    
%end

% --- LWP  ---
var='TGCLDLWP'; stash = ['002_' var file_suffix];   %Total grid-box cloud LWP
SW_calcs_James_King_CESM_Nov2022_load_input_data_FUNC 
%lwp_all_sky_levs = 1e3*dat; %convert to g/m2 from kg/m2
output.lwp_ic = dat(1:nT,:,:) ./ output.totCF;
output.lwp_ic(output.totCF<0.05) = NaN;



% Nd
% stash/dir 01241 - CDNC * layer weight
%UM_run = 'u-ck596'; %control run

%CDNC
% stash = '1241'; var='UM_m01s01i241_vn1105';
% SW_calcs_James_Weber_Jan2022_load_input_data_FUNC
% Nd_times_wt = dat; %cm^-3

var='NUMLIQ'; stash = ['002_' var file_suffix];   %Grid-box average liquid number (per kg)
SW_calcs_James_King_CESM_Nov2022_load_input_data_FUNC
Nd_4D = dat(1:nT,:,:,:)/1e6; %convert to cm^-3
lev=output.plevs;

%Temperature
%var='T'; stash = ['002_' var file_suffix];   %Grid-box average liquid number (per kg)
%SW_calcs_James_King_CESM_Nov2022_load_input_data_FUNC

%air density
var='RHO_CLUBB'; stash = ['002_' var file_suffix];   %Grid-box average liquid number (per kg)
SW_calcs_James_King_CESM_Nov2022_load_input_data_FUNC
%RHO_CLUBB is on 33 levels that are different to the levels that Nd is on.
%I'm presuming here that the output is indeed on pressure levels. And then
%I'm interpolating each column of RHO_CLUBB onto the lev for Nd. Not sure
%if the "atmosphere_hybrid_sigma_pressure_coordinate" in hPa is indeed a
%pressure level, or if it needs converting to actual pressure using a
%formula.
dat=permute(dat(1:nT,:,:,:),[2 3 4 1]); %permute to put the height column at the end for interp1 to work on a matrix
ilev = output.pilevs;
ilev_bnds = output.pilevs_bnds;
rhoa = interp1(ilev,dat,lev);
rhoa = permute(rhoa,[4 1 2 3 ]);
Nd_4D = Nd_4D .* rhoa; %convert to #/cm3 from #/kg

% Liquid cloud frequency for calculation of in-cloud CDNC
var='FREQL'; stash = ['002_' var file_suffix];
SW_calcs_James_King_CESM_Nov2022_load_input_data_FUNC

%Divide by the liquid cloud fraction to get the in-cloud Nd
output.Nd_ic = Nd_4D(1:nT,:,:,:)./ dat(1:nT,:,:,:); %
output.Nd_ic(dat(1:nT,:,:,:)<0.05) = NaN;


% Grid-box mean LWC
var='LWC'; stash = ['002_' var file_suffix];
SW_calcs_James_King_CESM_Nov2022_load_input_data_FUNC

%Verticaly integrate in-cloud CDNC weighted by LWC
%plevs = unique(output.plevs_bnds(:));
%dz = diff(plevs);
%Don't bother with dz for now actually.
[output.Nd,N,std_dev] = meanNoNan(output.Nd_ic,2,'',1,1,dat(1:nT,:,:,:));


output.lat2d = output.gcm_Plat2D_UM;
output.lat2d_edges = output.gcm_Plat2D_edges_UM;
output.lon2d = output.gcm_Plon2D_UM;
output.lon2d_edges = output.gcm_Plon2D_edges_UM;

% --- Aerosol optical depth


% --- Absorbing aerosol optical depth


% --- Dust optical depth



end

