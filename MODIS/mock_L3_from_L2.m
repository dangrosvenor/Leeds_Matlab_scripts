try
    %think here we need to specify the edges of the box that we want, since
    %finds points within the values specified by LATS and LONS


    %this groups data in one degree lat/lon areas
    %see save_mockL3_vars.m for a routine to save these fields
    %need one for loading?

    %calculate Plat_L2, Plon_L2, etc.
    %filtering_data_L2; - now done in open_L2_MODIS_file_01

    %choose to make the arrays where each data individual pixel within each lat/lon point is
    %stored - probably not needed (and uses too much memory) for processing many files
    i_make_all=0;
    
    ireturn=0;

    itime=1;
    NTL2 = 1; %size of the time dimension

    thresh_Nd_per_error=100;

    if exist('override_mockL3_options') & override_mockL3_options==1
        clear override_mockL3_options
    else
%        box_type = 'NxN pixel square';
%        Npix=5; %makes squares of Npix*Npix

        box_type = 'regular lat lon grid';
        LAT_step=1;
        LON_step=1;

    end

    if ~exist('override_mockL3_options2') | override_mockL3_options2==0

        %choose the lat/lon limits for the region we want to process
        LAT_min = floor(minALL(Plat_L2));
        LAT_max=  ceil(maxALL(Plat_L2));

        LON_min = floor(minALL(Plon_L2));
        LON_max = ceil(maxALL(Plon_L2));

        LAT_min=67; LAT_max=75; LON_min=-180; LON_max=-135;  %MPACE Arctic
        %LAT_min=-85; LAT_max=-40; LON_min=-180; LON_max=180; %Antarctica
        %LAT_min=-89; LAT_max=-40; LON_min=-180; LON_max=0; %Antarctica
        LAT_min=40; LAT_max=89; LON_min=-90; LON_max=90; %Arctic summer 2004 files
        LAT_min=40; LAT_max=90; LON_min=-180; LON_max=180; %Arctic summer 2007 files
        %LAT_min=0; LAT_max=90; LON_min=-157; LON_max=-149; %Arctic summer 2004 files - for max SZA - just for the region
        %that was downloaded
        LAT_min = 72; LAT_max=75; LON_min=-3; LON_max=48; %Arctic summer
        %LAT_min=-23; LAT_max=-16; LON_min=-90; LON_max=-80; %VOCALS Zhang & Plantnick test region MOD06_L2.A2005092.1600.051.2010295045517.hdf

        %LAT_min=-90; LAT_max=90; LON_min=-180; LON_max=180; %Global


    end




    switch box_type
        case 'Selected NxN pixel region';
            %just select one lat and lon
            LATS=[lat_new; lat_new];
            LONS=[lon_new; lon_new];
            mockL3_box_type_str=['Selected ' num2str(Npix) 'x' num2str(Npix) ' pixel region '];

            i_LATS=1; %just will do one loop
            i_LONS=1;

        case 'NxN pixel square';
            a=[1:Npix:size(Plat_L2,1)];
            b=[1:Npix:size(Plat_L2,2)];
            LATS=1:length(a);
            LONS=1:length(b);
            mockL3_box_type_str=[num2str(Npix) 'x' num2str(Npix) ' pixel square '];
        otherwise
            mockL3_box_type_str=[num2str(LAT_step) 'x' num2str(LAT_step) ' degree square'];
            LATS=[LAT_min:LAT_step:LAT_max];
            LONS=[LON_min:LON_step:LON_max];
    end



    switch box_type
        case 'Selected NxN pixel region'
            %For this case choose the pixels now as we will only have one value
            %to calculate and can exit if we don't have this region in the
            %swath

            dist_tol = 3; % km tolerance in the distance to the lat lon

            %here we want an Npix x Npix square centred on a specific
            %location (lat_new, lon_new)

            %LATS and LONS are just one location now.
            %indices (referencing the large array) of the points in the required lat/lon point

            %find the approx region
            iapprox=find(Plat_L2>=LATS(1)-2 & Plat_L2<LATS(1)+2 & Plon_L2>=LONS(1)-2 & Plon_L2<LONS(1)+2);

            if length(iapprox)>0

                [dist_NN] = distlatlon(Plat_L2(iapprox),Plon_L2(iapprox),lat_new,lon_new);
                [dist_min,idist] = min(dist_NN);

                if isnan(dist_min) | dist_min > dist_tol
                    ireturn=1; %exit this subroutine
                else
                    %choose the centre of the selected pixel - better
                    %to have Npix an odd number

                    [xpos,ypos] = ind2sub(size(Plat_L2),iapprox(idist));
                    di = (Npix-1)/2;
                    ix=[xpos-di:xpos+di];
                    iy=[ypos-di:ypos+di];
                    [ix2,iy2]=meshgrid(ix,iy);
                    %linear indices
                    sLAT_LON = size(Plat_L2);

                    if min(ix)>=1 & max(ix)<=sLAT_LON(1) & min(iy)>=1 & max(iy)<=sLAT_LON(2)
                        ipos_L2 = sub2ind(size(Plat_L2),ix2(:),iy2(:));
                    else
                        ireturn=1; %indices for the square are out of the range of the swath
                    end


                end

            else
                ireturn=1;

            end


            %  -- Should allow for instances when the requested point lies at the edge of the swath since the
            % requested Npix x Npix box may cross into the next swath. Could fix
            % this by loading in two swaths at once and joining them together.


        otherwise

            %% only process the lat/lon values that are actually contained within the
            %swath to save time

            [i1,i2]=findheight_nearest(LATS,floor(minALL(Plat_L2)),ceil(maxALL(Plat_L2)),'bounded');
            if (i1==0 & i2==0) | (i1>1e19 & i2>1e19) %if this is the case then the whole range is out of our LATS bound
                i_LATS=[]; %if have no values in range then don't process
            else  %i1=0 when some of the requested range is outside of LATS(1), but some is inside. Simlar for i2=1e20
                if i1==0
                    i1=1;
                end
                if i2>1e19;
                    i2=length(LATS);
                end

                i_LATS = i1:min(i2,length(LATS)-1); %need to take one off as the full loop is 1:length(LATS)-1

            end

            [i1,i2]=findheight_nearest(LONS,floor(minALL(Plon_L2)),ceil(maxALL(Plon_L2)),'bounded');
            if (i1==0 & i2==0) | (i1>1e19 & i2>1e19) %if this is the case the whole range is out of our LOiNS bound
                i_LONS=[]; %if have no values in range then don't process
            else
                if i1==0
                    i1=1;
                end
                if i2>1e19;
                    i2=length(LONS);
                end
                i_LONS = i1:min(i2,length(LONS)-1); %need to take one off as the full loop is 1:length(LONS)-1
            end
            %ipos_L2=find(Plat_L2>=LATS(ilatL3L2) & Plat_L2<LATS(ilatL3L2+1) & Plon_L2>=LONS(ilonL3L2) & Plon_L2<LONS(ilonL3L2+1));



    end


    %in case we just want to know the size of the arrays to be stored
    %(e.g. as for MODIS_process_multiple_L2_files
    if exist('ijust_size_check') & ijust_size_check==1
        clear ijust_size_check
        return
    end



%% set up the screening required - points to be removed. Note the NOT sign in
%the find brackets

                                
%                                ihtot = find( ~ ( cfL2>=thresh_CF ) ); thresh_str_mock_L3=['CF.GTE.' num2str(thresh_CF)];
%                                ihtot = find( ~ ( percent_error_Nd<thresh_Nd_per_error ) ); thresh_str_mock_L3=['Nd Perror.LT.' num2str(thresh_Nd_per_error)];
%                                ihtot = find( ~ ( re_un<thresh_Reff_per_error ) ); thresh_str_mock_L3=['Reff Perror.LT.' num2str(thresh_Reff_per_error)];                                
%                                ihtot = find( ~ ( re_un_abs<thresh_Reff_abs_error ) ); thresh_str_mock_L3=['Reff abs error.LT.' num2str(thresh_Reff_abs_error)];                                                                
%                                ihtot = find( ~ ( percent_error_Nd<thresh_Nd_per_error & N>400) ); thresh_str_mock_L3=['Perror.LT.' num2str(thresh_Nd_per_error) '.AND.N.GT.400'];
%                                ihtot = find( ~ ( re_un<thresh_Reff_per_error & (phase_flag==2) ) ); thresh_str_mock_L3=['Reff % error.LT.' num2str(thresh_Reff_per_error) ' AND liquid phase.'];  %                                 ihtot = find( ~ ( percent_error_Nd<thresh_Nd_per_error & re<thresh_Reff) ); thresh_str_mock_L3=['Perror.LT.' num2str(thresh_Nd_per_error) 'AND Re.LT.' num2str(thresh_Reff)];
%                                ihtot = find( ~ ( percent_error_Nd<thresh_Nd_per_error & (phase_flag==2) ) ); thresh_str_mock_L3=['Nd % error.LT.' num2str(thresh_Nd_per_error) ' AND liquid phase.'];  
%                                ihtot = find( ~ ( re_un_abs<thresh_Reff_abs_error & (phase_flag==3) ) ); thresh_str_mock_L3=['Reff abs error.LT.' num2str(thresh_Reff_abs_error) ' AND ice phase.'];                                  
%                                ihtot = find( ~ ( phase_flag==2 ) ); thresh_str_mock_L3=[' liquid phase.']; 
%                                ihtot = find( ~ ( phase_flag==3 & tau_bounds==0 )  ); thresh_str_mock_L3=[' ICE phase.AND.Tau within bounds'];                                 
%                                 ihtot = find( ~ ( phase_flag==2 & tau_bounds==0) ); thresh_str_mock_L3=[' liquid phase.AND.Tau within bounds'];

%                                 ihtot = find( ~ ( phase_flag==2 & tau_bounds==0 & surface_flag==0) ); thresh_str_mock_L3=[' liquid phase.AND.Tau within bounds, ocean only'];
%                                 ihtot = find( ~ ( CTT>thresh_CTT ) ); thresh_str_mock_L3=[' CTT.GT.' num2str(thresh_CTT)]; 
%                                 ihtot = find( ~ ( percent_error_Nd<thresh_Nd_per_error & phase_flag==2 & tau_bounds==0) ); thresh_str_mock_L3=['Perror.LT.' num2str(thresh_Nd_per_error) ' ,liquid phase.AND.Tau within bounds'];                                 
%                                ihtot = find( ~ ( phase_flag==2 & squeeze(mask_1km(1,:,:))==1) ); thresh_str_mock_L3=[' liquid phase, cloud mask determined'];   
%                                ihtot = find( ~ ( phase_flag==2 & squeeze(mask_1km(1,:,:))==1 & squeeze(mask_1km(4,:,:)) == 1 )  ); thresh_str_mock_L3=[' liquid phase, cloud mask determined, no sunglint']; 
                                ihtot = find( ~ ( phase_retreival_outcome==1 & phase_flag==2 & squeeze(mask_1km(4,:,:)) == 1  )  ); thresh_str_mock_L3=[' liquid phase, successful phase outcome, no sunglint'];                                   
%                                ihtot = find( ~ ( phase_retreival_outcome==1 & phase_flag==2 & squeeze(mask_1km(4,:,:)) == 1  & squeeze(qapp_1km(16,:,:)) == 2 )  ); thresh_str_mock_L3=[' liquid phase, successful phase outcome, no sunglint, single liquid layer'];
                                
                                
                                %mask_1km(1,:,:)==0 means that the cloud mask was undetermined
                                 % mask_1km(4,:,:) - sunglint, 0 = yes (i.e. bad data), 1 = no (good data)
                                 
%                                   ihtot=[]; thresh_str_mock_L3='xxx'; %set this for no screening
 %indices for CTT in order to calculate CTT for ice, liquid and undetermined phases
 ihtot_all = find( ~ ( phase_retreival_outcome==1 & squeeze(mask_1km(4,:,:)) == 1  )  ); thresh_str_mock_L3_2=[' successful phase outcome, no sunglint'];                                                                   

%% Set up arrays, etc. 
                                   
%when looking at re for the different wavelengths note that have made NaN values when the cloud retrieval
%was not confident - i.e. when made re values NaN also made the re_diff values
%NaN in open_L2_MODIS_file_01.m (or open_L2_Joint_MODIS_file_02)
%But, it seems that for 1.6 um we get many low Re values (down to zero),
%whereas 2.1 um Re values tend to not go below 5um
re_min=0.5;
re_diff_max=15;
                                   
re_21 = re;
re_21(ihtot)=NaN;
re_21(re_21<re_min)=NaN;
re_16 = re_21 + squeeze(re_diff(:,:,1));
re_37 = re_21 + squeeze(re_diff(:,:,2));

%Note - now for any pixels where re21 is NaN the other wavelengths will be
%NaN since they are calculated from re21 + re_diff.
%E.g. if confidence screening makes re21 values NaN then the others will be
%NaN too. Although there will be instances when re37 etc will go over 20um
%(sincre re21 might be below 20um).


%calculate the droplet field and do the screening
Wflag='calc';
[N,H,W,k,Q,cw]=MODIS_N_H_func(tau,re_21*1e-6,Wflag,NaN,t_top2); %added CTT
N(ihtot)=NaN; %make these values NaN as they will then be removed from the average
W=W*1e3; %convert to g/m2
W(ihtot)=NaN;



[N_16,H_16,W_16,k,Q,cw]=MODIS_N_H_func(tau,re_16*1e-6,Wflag,NaN,t_top2); %added CTT
if ~exist('i_limit_16_37_to_20_um') | i_limit_16_37_to_20_um==0
    ihtot2 = find (abs(re_16-re_21)>re_diff_max | re_16<re_min);    
else
    ihtot2 = find (abs(re_16-re_21)>re_diff_max | re_16<re_min | re_16>20);
end
ihtot2=union(ihtot,ihtot2);
N_16(ihtot2)=NaN; %make these values NaN as they will then be removed from the average
W_16=W_16*1e3; %convert to g/m2
W_16(ihtot2)=NaN; 
re_16(ihtot2)=NaN;

[N_37,H_37,W_37,k,Q,cw]=MODIS_N_H_func(tau,re_37*1e-6,Wflag,NaN,t_top2); %added CTT
if ~exist('i_limit_16_37_to_20_um') | i_limit_16_37_to_20_um==0
    ihtot2 = find (abs(re_37-re_21)>re_diff_max | re_37<re_min);
else
    %Here will be limiting re16 and re37 to 20um to match what happens to
    %re21 when we apply 'very good' only water path confidence screening
    %(when combined with liquid clouds).
    ihtot2 = find (abs(re_37-re_21)>re_diff_max | re_37<re_min | re_37>20 );
end
ihtot2=union(ihtot,ihtot2);
N_37(ihtot2)=NaN; %make these values NaN as they will then be removed from the average
W_37=W_37*1e3; %convert to g/m2
W_37(ihtot2)=NaN; 
re_37(ihtot2)=NaN;

%Find pixels where any of the wavelengths are NaN and make those NaN for
%all wavelengths to give consistency between the three (all usign the same
%pixels)
if exist('i_make_wavelengths_consistent') & i_make_wavelengths_consistent==1
    inan_any = find( isnan(re_37)==1 | isnan(re_21)==1 | isnan(re_16)==1 );
    
    N_16(inan_any)=NaN; %make these values NaN as they will then be removed from the average
    W_16(inan_any)=NaN;
    re_16(inan_any)=NaN;    
    
    N_37(inan_any)=NaN; %make these values NaN as they will then be removed from the average
    W_37(inan_any)=NaN;
    re_37(inan_any)=NaN;
    
    N(inan_any)=NaN; %make these values NaN as they will then be removed from the average
    W(inan_any)=NaN;
    re(inan_any)=NaN;    
    
end

%choose the regions which *I think* correspond to the lat/lon grid made
%from the 5km grid. Note - these are now set in filtering_data_L2 using the
%correct swath sampling as given from the variable attributes
% row_L2_inds = [4:1349];
% col_L2_inds = [4:2029];

%note that the values here have been screened by ihtot already
Nd_L2_swath = N(row_L2_inds,col_L2_inds);
Nd_L2_swath_16 = N_16(row_L2_inds,col_L2_inds);
Nd_L2_swath_37 = N_37(row_L2_inds,col_L2_inds);
W_L2_swath = W(row_L2_inds,col_L2_inds);
W_L2_swath_16 = W_16(row_L2_inds,col_L2_inds);
W_L2_swath_37 = W_37(row_L2_inds,col_L2_inds);
cfL2_2 = cfL2(row_L2_inds,col_L2_inds);

Tau_L2_swath = tau(row_L2_inds,col_L2_inds);
Re_L2_swath = re_21(row_L2_inds,col_L2_inds);
Re_L2_swath_16 = re_16(row_L2_inds,col_L2_inds);
Re_L2_swath_37 = re_37(row_L2_inds,col_L2_inds);
phase_flag_2 = phase_flag(row_L2_inds,col_L2_inds); %N.B. phase_flag = squeeze(qapp_1km(1,:,:)); 
%      [qapq_1km,qapp_1km] = flagread_1km_Dan(qa1,1);
%      qapp_1km is a continuation of the 1km QA as documented in the MODIS QA
%      plan document. Starts at the 3rd byte ("Primary Cloud Retrieval Phase
%      Flag").
phase_retreival_outcome_2 = phase_retreival_outcome(row_L2_inds,col_L2_inds);
%      phase_retreival_outcome = squeeze(qapp_1km(2,:,:)); %Outcome of primary phase retrieval above. 0=not attempted
%      (possibly clear, although can actually have a phase for these - but it
%      should be ignored). 1=successful, is cloudy. All clear points have =0
%      (&phase_flag==1)
sunglint_flag = squeeze(mask_1km(4,row_L2_inds,col_L2_inds)); 

multi_layer_flag = squeeze(qapp_1km(16,row_L2_inds,col_L2_inds));

tau_bounds2 = tau_bounds(row_L2_inds,col_L2_inds);


%don't think we need to make these NaN for ihtot since only select data for
%inclusion below (using igood2) when Nd is not NaN

tau_useful = squeeze(qapq_1km(1,:,:)); tau_useful = tau_useful(row_L2_inds,col_L2_inds);
tau_confidence = squeeze(qapq_1km(2,:,:)); tau_confidence = tau_confidence(row_L2_inds,col_L2_inds);
mask_det = squeeze(mask_1km(1,:,:)); mask_det = mask_det(row_L2_inds,col_L2_inds);
cloudy = squeeze(mask_1km(2,:,:)); cloudy = cloudy(row_L2_inds,col_L2_inds);
Tau_un_L2_swath = tau_un;           Tau_un_L2_swath(ihtot)=NaN; Tau_un_L2_swath = Tau_un_L2_swath(row_L2_inds,col_L2_inds);
Re_un_L2_swath =  re_un;            Re_un_L2_swath(ihtot)=NaN;  Re_un_L2_swath =  Re_un_L2_swath(row_L2_inds,col_L2_inds);
Nd_un_L2_swath =  percent_error_Nd; Nd_un_L2_swath(ihtot)=NaN;  Nd_un_L2_swath =  Nd_un_L2_swath(row_L2_inds,col_L2_inds);
%CTT_L2_swath = t_top2;           CTT_L2_swath(ihtot)=NaN; CTT_L2_swath = CTT_L2_swath(row_L2_inds,col_L2_inds);
%Have removed the screening of CTT for liquid only points here since we
%want to look at CTT in ice regions too.
CTT_L2_swath = t_top2;  CTT_L2_swath = CTT_L2_swath(row_L2_inds,col_L2_inds);
CTP_L2_swath = p_top2;           CTP_L2_swath(ihtot)=NaN; CTP_L2_swath = CTP_L2_swath(row_L2_inds,col_L2_inds);
%solarZA_L2_swath = solarZA_L2;    solarZA_L2_swath(ihtot)=NaN;  solarZA_L2_swath=solarZA_L2_swath(row_L2_inds,col_L2_inds);
%solarAz_L2_swath = solarAz_L2;    solarAz_L2_swath(ihtot)=NaN;  solarAz_L2_swath=solarAz_L2_swath(row_L2_inds,col_L2_inds);
%sensorZA_L2_swath = sensorZA_L2;    sensorZA_L2_swath(ihtot)=NaN;  sensorZA_L2_swath=sensorZA_L2_swath(row_L2_inds,col_L2_inds);
%sensorAz_L2_swath = sensorAz_L2;    sensorAz_L2_swath(ihtot)=NaN;  sensorAz_L2_swath=sensorAz_L2_swath(row_L2_inds,col_L2_inds);
%removed the setting of these to NaN as we want them for all points, not
%just the ones with liquid cloud
solarZA_L2_swath = solarZA_L2;      solarZA_L2_swath=solarZA_L2_swath(row_L2_inds,col_L2_inds);
solarAz_L2_swath = solarAz_L2;      solarAz_L2_swath=solarAz_L2_swath(row_L2_inds,col_L2_inds);
sensorZA_L2_swath = sensorZA_L2;    sensorZA_L2_swath=sensorZA_L2_swath(row_L2_inds,col_L2_inds);
sensorAz_L2_swath = sensorAz_L2;    sensorAz_L2_swath=sensorAz_L2_swath(row_L2_inds,col_L2_inds);


%the bins
Nd_bins=[0:20:5000];
W_bins=[0:20:6000];
Re_bins=[0:0.4:100];
Tau_bins=[0:0.4:100];
Nd_un_bins=[0:1.2:300];
Tau_un_bins=[0:1.2:300];
Re_un_bins=[0:1.2:300];

Tau2D_bins = Tau_bins;
Re2D_bins = Re_bins;

Tau2D_bins = [0:1.6:100];
Re2D_bins = [0:1.6:100];

Nd_2D_bins = [0:5:2000];
CTP_2D_bins = [120:25:1020];


%create arrays to start with to be more memory efficient and make all NaN
meanNd_mockL3 = NaN*ones(length(LATS)-1,length(LONS)-1);
meanNd_mockL3_16 = NaN*ones(length(LATS)-1,length(LONS)-1);
meanNd_mockL3_37 = NaN*ones(length(LATS)-1,length(LONS)-1);
meanTau_mockL3 = NaN*ones(length(LATS)-1,length(LONS)-1);
meanRe_mockL3 = NaN*ones(length(LATS)-1,length(LONS)-1);
meanRe_mockL3_16 = NaN*ones(length(LATS)-1,length(LONS)-1);
meanRe_mockL3_37 = NaN*ones(length(LATS)-1,length(LONS)-1);
meanRe_log_mockL3 = NaN*ones(length(LATS)-1,length(LONS)-1);
meanRe_log_mockL3_16 = NaN*ones(length(LATS)-1,length(LONS)-1);
meanRe_log_mockL3_37 = NaN*ones(length(LATS)-1,length(LONS)-1);
meanNd_un_mockL3 = NaN*ones(length(LATS)-1,length(LONS)-1);
meanNd_un_combined_mockL3 = NaN*ones(length(LATS)-1,length(LONS)-1);
meanTau_un_mockL3 = NaN*ones(length(LATS)-1,length(LONS)-1);
meanTau_log_mockL3 = NaN*ones(length(LATS)-1,length(LONS)-1);
meanTau_logun_mockL3 = NaN*ones(length(LATS)-1,length(LONS)-1);
meanRe_un_mockL3 = NaN*ones(length(LATS)-1,length(LONS)-1);
meanRe_logun_mockL3 = NaN*ones(length(LATS)-1,length(LONS)-1);
meanCTT_mockL3 = NaN*ones(length(LATS)-1,length(LONS)-1);
stdCTT_mockL3 = NaN*ones(length(LATS)-1,length(LONS)-1);
meanCTT_ALL_mockL3 = NaN*ones(length(LATS)-1,length(LONS)-1); %new values that give CTT for all pixels (no liquid screening)
stdCTT_ALL_mockL3 = NaN*ones(length(LATS)-1,length(LONS)-1);
meanCTT_CF_mockL3 = NaN*ones(length(LATS)-1,length(LONS)-1); %new values that give CTT for all pixels used for CF_liq
stdCTT_CF_mockL3 = NaN*ones(length(LATS)-1,length(LONS)-1);

minCTT_mockL3 = NaN*ones(length(LATS)-1,length(LONS)-1);
maxCTT_mockL3 = NaN*ones(length(LATS)-1,length(LONS)-1);
meanCTT_ice_mockL3 = NaN*ones(length(LATS)-1,length(LONS)-1);
minCTT_ice_mockL3 = NaN*ones(length(LATS)-1,length(LONS)-1);
maxCTT_ice_mockL3 = NaN*ones(length(LATS)-1,length(LONS)-1);
meanSolarZA_mockL3 = NaN*ones(length(LATS)-1,length(LONS)-1);
meanSolarAz_mockL3 = NaN*ones(length(LATS)-1,length(LONS)-1);
meanSensorZA_mockL3 = NaN*ones(length(LATS)-1,length(LONS)-1);
meanSensorAz_mockL3 = NaN*ones(length(LATS)-1,length(LONS)-1);
meanW_mockL3 = NaN*ones(length(LATS)-1,length(LONS)-1);
stdW_mockL3 = NaN*ones(length(LATS)-1,length(LONS)-1);
stdTau_mockL3 = NaN*ones(length(LATS)-1,length(LONS)-1);
stdRe_mockL3 = NaN*ones(length(LATS)-1,length(LONS)-1);
stdRe_mockL3_16 = NaN*ones(length(LATS)-1,length(LONS)-1);
stdRe_mockL3_37 = NaN*ones(length(LATS)-1,length(LONS)-1);
minRe_mockL3 = NaN*ones(length(LATS)-1,length(LONS)-1);
minRe_mockL3_16 = NaN*ones(length(LATS)-1,length(LONS)-1);
minRe_mockL3_37 = NaN*ones(length(LATS)-1,length(LONS)-1);
maxRe_mockL3 = NaN*ones(length(LATS)-1,length(LONS)-1);
maxRe_mockL3_16 = NaN*ones(length(LATS)-1,length(LONS)-1);
maxRe_mockL3_37 = NaN*ones(length(LATS)-1,length(LONS)-1);
minTau_mockL3 = NaN*ones(length(LATS)-1,length(LONS)-1);
maxTau_mockL3 = NaN*ones(length(LATS)-1,length(LONS)-1);


stdNd_mockL3 = NaN*ones(length(LATS)-1,length(LONS)-1);
stdNd_mockL3_16 = NaN*ones(length(LATS)-1,length(LONS)-1);
stdNd_mockL3_37 = NaN*ones(length(LATS)-1,length(LONS)-1);
meanlogW_mockL3 = NaN*ones(length(LATS)-1,length(LONS)-1);

LAT_mockL3 = NaN*ones(length(LATS)-1,length(LONS)-1);
LON_mockL3 = NaN*ones(length(LATS)-1,length(LONS)-1);
%edge array needs to be one size larger
LAT_mockL3_edge = NaN*ones(length(LATS),length(LONS));
LON_mockL3_edge = NaN*ones(length(LATS),length(LONS));

CF_mockL3 = NaN*ones(length(LATS)-1,length(LONS)-1,NTL2); %Cloud fraction
Np_mockL3 = NaN*ones(length(LATS)-1,length(LONS)-1,NTL2); %number of datapoints used
Nptot_mockL3 = NaN*ones(length(LATS)-1,length(LONS)-1,NTL2); %total no. possible datapoints
Npmask_mockL3 = NaN*ones(length(LATS)-1,length(LONS)-1,NTL2); %number of datapoints used

Nice_mockL3 = NaN*ones(length(LATS)-1,length(LONS)-1,NTL2); %number of pixels identified as ice
Nice_mockL3_new = NaN*ones(length(LATS)-1,length(LONS)-1,NTL2); %
Nliq_mockL3 = NaN*ones(length(LATS)-1,length(LONS)-1,NTL2); %
Nliq_mockL3_new = NaN*ones(length(LATS)-1,length(LONS)-1,NTL2); %
N_clear_new = NaN*ones(length(LATS)-1,length(LONS)-1,NTL2); %

Nundet_mockL3 = NaN*ones(length(LATS)-1,length(LONS)-1,NTL2); %
Nundet_mockL3_new = NaN*ones(length(LATS)-1,length(LONS)-1,NTL2); %

N_no_opt_or_clear_mockL3 = NaN*ones(length(LATS)-1,length(LONS)-1,NTL2); %
Nreject_mockL3 = NaN*ones(length(LATS)-1,length(LONS)-1,NTL2); %
N_cloudy_confident_mockL3 = NaN*ones(length(LATS)-1,length(LONS)-1,NTL2); %
N_cloudy_probable_mockL3 = NaN*ones(length(LATS)-1,length(LONS)-1,NTL2); %
N_cloudy_no_opt_mockL3 = NaN*ones(length(LATS)-1,length(LONS)-1,NTL2); %

N_single_liquid_layer_mockL3 = NaN*ones(length(LATS)-1,length(LONS)-1,NTL2); %

NdPDF_mockL3 = NaN*ones(length(LATS)-1,length(LONS)-1,length(Nd_bins)-1,NTL2); %1D PDF of Nd
NdPDF_mockL3_16 = NaN*ones(length(LATS)-1,length(LONS)-1,length(Nd_bins)-1,NTL2); %1D PDF of Nd
NdPDF_mockL3_37 = NaN*ones(length(LATS)-1,length(LONS)-1,length(Nd_bins)-1,NTL2); %1D PDF of Nd
WPDF_mockL3 = NaN*ones(length(LATS)-1,length(LONS)-1,length(W_bins)-1,NTL2); %1D PDF of Nd
WPDF_mockL3_16 = NaN*ones(length(LATS)-1,length(LONS)-1,length(W_bins)-1,NTL2); %1D PDF of Nd
WPDF_mockL3_37 = NaN*ones(length(LATS)-1,length(LONS)-1,length(W_bins)-1,NTL2); %1D PDF of Nd
RePDF_mockL3 = NaN*ones(length(LATS)-1,length(LONS)-1,length(Re_bins)-1,NTL2); %1D PDF of Nd
RePDF_mockL3_16 = NaN*ones(length(LATS)-1,length(LONS)-1,length(Re_bins)-1,NTL2); %1D PDF of Nd
RePDF_mockL3_37 = NaN*ones(length(LATS)-1,length(LONS)-1,length(Re_bins)-1,NTL2); %1D PDF of Nd
TauPDF_mockL3 = NaN*ones(length(LATS)-1,length(LONS)-1,length(Tau_bins)-1,NTL2); %1D PDF of Nd

Nd_un_PDF_mockL3 = NaN*ones(length(LATS)-1,length(LONS)-1,length(Nd_un_bins)-1,NTL2); %1D PDF of Nd
Re_un_PDF_mockL3 = NaN*ones(length(LATS)-1,length(LONS)-1,length(Re_un_bins)-1,NTL2); %1D PDF of Nd
Tau_un_PDF_mockL3 = NaN*ones(length(LATS)-1,length(LONS)-1,length(Tau_un_bins)-1,NTL2); %1D PDF of Nd

Nd_vs_CTP_mockL3 = NaN*ones(length(LATS)-1,length(LONS)-1,length(CTP_2D_bins)-1,NTL2); %1D PDF of Nd
NP_vs_CTP_mockL3 = NaN*ones(length(LATS)-1,length(LONS)-1,length(CTP_2D_bins)-1,NTL2); %1D PDF of Nd

Tau_Re_2D_PDF = NaN*ones(length(LATS)-1,length(LONS)-1,length(Tau2D_bins)-1,length(Re2D_bins)-1,NTL2); %2D PDF of Tau and Reff

%estimate the max possible number of pixels within each cell (approx
%110*110km at 1km pixel size)
Ncell_max=120*120; 

if i_make_all==1

    %create an array to store all of the exact Nd values in
    Nd_all = NaN*ones(length(LATS)-1,length(LONS)-1,Ncell_max);
    Tau_all = NaN*ones(length(LATS)-1,length(LONS)-1,Ncell_max);
    Re_all = NaN*ones(length(LATS)-1,length(LONS)-1,Ncell_max);
    Tau_un_all = NaN*ones(length(LATS)-1,length(LONS)-1,Ncell_max);
    Re_un_all = NaN*ones(length(LATS)-1,length(LONS)-1,Ncell_max);
    Nd_un_all = NaN*ones(length(LATS)-1,length(LONS)-1,Ncell_max);
    CTT_all = NaN*ones(length(LATS)-1,length(LONS)-1,Ncell_max);
    SolarZA_all = NaN*ones(length(LATS)-1,length(LONS)-1,Ncell_max);

end

% if we want to exit - doing so after we've made the NaN arrays, so that
% NaN values will be returned rather than the last value
if ireturn==1
    return
end


%% Main loop

    for ilonL3L2=i_LONS %swapped these so that lat index (first one) is varying quickest
%        ilonL3L2
        for ilatL3L2=i_LATS   %in order to follow the way it is ordered in memory


        
        switch box_type
            case 'Selected NxN pixel region'
                %ipos_L2 calculated earlier
                
                LAT_mockL3(ilatL3L2,ilonL3L2,itime) = mean(Plat_L2(ipos_L2),1);
                LON_mockL3(ilatL3L2,ilonL3L2,itime) = mean(Plon_L2(ipos_L2),1);
                %want the edges of the averaging boxes for plotting purposes -
                %should be the first pixel - or could do an average along the box
                %edge?
                %        LAT_mockL3_edge(ilatL3L2,ilonL3L2,itime) = mean(Plat_L2(ix,iy(1)));
                %        LAT_mockL3_edge(ilatL3L2,ilonL3L2,itime) = mean(Plon_L2(ix(1),iy));

                LAT_mockL3_edge(ilatL3L2,ilonL3L2,itime) = Plat_L2(ix(1),iy(1));
                LON_mockL3_edge(ilatL3L2,ilonL3L2,itime) = Plon_L2(ix(1),iy(1));
                
                
            case 'NxN pixel square'
                %here we are dealing with the reduced grid of the same size
                %as Plat_L2 (after removing some of the swath pixels)
                xpos=(ilatL3L2-1)*Npix+1;
                ypos=(ilonL3L2-1)*Npix+1;                
                ix=[xpos:xpos+Npix-1];
                iy=[ypos:ypos+Npix-1];
                [ix2,iy2]=meshgrid(ix,iy);
                %linear indices
                ipos_L2 = sub2ind(size(Plat_L2),ix2(:),iy2(:));
                
                LAT_mockL3(ilatL3L2,ilonL3L2,itime) = mean(Plat_L2(ipos_L2),1);
                LON_mockL3(ilatL3L2,ilonL3L2,itime) = mean(Plon_L2(ipos_L2),1);
                %want the edges of the averaging boxes for plotting purposes -
                %should be the first pixel - or could do an average along the box
                %edge?
                %        LAT_mockL3_edge(ilatL3L2,ilonL3L2,itime) = mean(Plat_L2(ix,iy(1)));
                %        LAT_mockL3_edge(ilatL3L2,ilonL3L2,itime) = mean(Plon_L2(ix(1),iy));

                LAT_mockL3_edge(ilatL3L2,ilonL3L2,itime) = Plat_L2(ix(1),iy(1));
                LON_mockL3_edge(ilatL3L2,ilonL3L2,itime) = Plon_L2(ix(1),iy(1));
                
                
                
        
            case 'regular lat lon grid'

                %indices (referencing the large array) of the points in the required lat/lon point
                ipos_L2=find(Plat_L2>=LATS(ilatL3L2) & Plat_L2<LATS(ilatL3L2+1) & Plon_L2>=LONS(ilonL3L2) & Plon_L2<LONS(ilonL3L2+1));

        end
        
         
        
        
        %indices of points when we have good values of Nd (non-screened values)
        %within the small array 
        %this will be the number of cloudy points in our box (not total
        %pixels)
        igood=find(isnan(Nd_L2_swath(ipos_L2))==0);
        %indices of those points in the original large array
        igood2=ipos_L2(igood);
        
       
        
        %for mask_det 1 means was determined, 0 undetermined        
        Npmask_mockL3(ilatL3L2,ilonL3L2) = sum(mask_det(ipos_L2) );
        

        %number of good points - those within the required lat/lon and
        %which haven't been screened by ihtot
        Np_mockL3(ilatL3L2,ilonL3L2) = length(igood2);
        %total number of pixels in the lat lon square before screening
        Nptot_mockL3(ilatL3L2,ilonL3L2) = length(ipos_L2);
        
        
         %mean cloud fraction with the lat/lon cell. Using ipos_L2 here
            %since we want the cloud fraction within all possible pixels
            CF_mockL3(ilatL3L2,ilonL3L2,itime) = meanNoNan(cfL2_2(ipos_L2),1);
            %similarly for CTT - do parameters that include all points
            %(ice and liquid) in case we want to look at low cloud
            %fractions
            [meanCTT_ALL_mockL3(ilatL3L2,ilonL3L2,itime),Npoints_temp,stdCTT_ALL_mockL3(ilatL3L2,ilonL3L2,itime)] = meanNoNan(CTT_L2_swath(ipos_L2),1);
            ice=find(phase_flag_2(ipos_L2)==3);
            ice_new = find(phase_flag_2(ipos_L2)==3  &  phase_retreival_outcome_2(ipos_L2)==1 &  sunglint_flag(ipos_L2)==1);   
            
            undet=find(phase_flag_2(ipos_L2)==4);
            undet_new=find(phase_flag_2(ipos_L2)==4  &  phase_retreival_outcome_2(ipos_L2)==1 &  sunglint_flag(ipos_L2)==1);
            
            liq=find(phase_flag_2(ipos_L2)==2);
            liq_new = find( phase_flag_2(ipos_L2)==2 &  phase_retreival_outcome_2(ipos_L2)==1 &  sunglint_flag(ipos_L2)==1);
            
            if length(liq_new)>0
                [meanCTT_CF_mockL3(ilatL3L2,ilonL3L2,itime),Npoints_temp,stdCTT_CF_mockL3(ilatL3L2,ilonL3L2,itime)] = meanNoNan(CTT_L2_swath(liq_new),1);
            end
            
            no_opt=find(phase_flag_2(ipos_L2)==0); %when there was no optical retrieval, although
            %there could still be cloudy (or clear)
            
            
            % Make a cloud fraction from the mask for probably cloudy and
            % confident cloudy. Would need to add the two to get all cloud
            % considered as probable or confident.
            %mask_5km(1,:,:) - Cloud Mask Status Flag, 0=undetermined, 1=determined
            %mask_5km(2,:,:) - Cloud Mask Cloudiness Flag, 0=confident cloudy,
            %    1=probably cloudy, 2=probably clear, 3=confident clear
            %See http://modis-atmos.gsfc.nasa.gov/_docs/QA_Plan_2011_01_26.pdf
            cloud_confident=find(mask_det(ipos_L2)==1 & cloudy(ipos_L2)==0);
            cloud_probable=find(mask_det(ipos_L2)==1 & cloudy(ipos_L2)==1);                          
            
            cloudy_no_opt=find(mask_det(ipos_L2)==1 & cloudy(ipos_L2)==0 & phase_flag_2(ipos_L2)==0); %pixels where there
            %is a confident cloud retrieval, but no optical property
            %determination (no phase retrieal attempted). May be useful to
            %give an idea of general confidence in a grid square
            
            single_liquid_layer = find(multi_layer_flag(ipos_L2))==0;
            
            %all clear pixles have phase_retreival_outcome_2==0
            clear_new = find( phase_flag_2(ipos_L2)==1 &  phase_retreival_outcome_2(ipos_L2)==0 & sunglint_flag(ipos_L2)==1);
            
            Nice_mockL3(ilatL3L2,ilonL3L2,itime) = length(ice);
            Nice_mockL3_new(ilatL3L2,ilonL3L2,itime) = length(ice_new);
             
            Nundet_mockL3(ilatL3L2,ilonL3L2,itime) = length(undet);
            Nundet_mockL3_new(ilatL3L2,ilonL3L2,itime) = length(undet_new);
                        
            Nliq_mockL3(ilatL3L2,ilonL3L2,itime) = length(liq);
            Nliq_mockL3_new(ilatL3L2,ilonL3L2,itime) = length(liq_new);
            
            N_clear_new(ilatL3L2,ilonL3L2,itime) = length(clear_new);
            
            N_no_opt_or_clear_mockL3(ilatL3L2,ilonL3L2,itime) = length(no_opt);  
            N_cloudy_confident_mockL3(ilatL3L2,ilonL3L2,itime) = length(cloud_confident); 
            N_cloudy_probable_mockL3(ilatL3L2,ilonL3L2,itime) = length(cloud_probable);             
            N_cloudy_no_opt_mockL3(ilatL3L2,ilonL3L2,itime) = length(cloudy_no_opt);  
            
            N_single_liquid_layer_mockL3(ilatL3L2,ilonL3L2,itime) = length(single_liquid_layer);              
            
            
            %using this: Cloud_Fraction_Liquid.timeseries3(:,:,iswath) =
            %  Nliq_mockL3 ./ Npmask_mockL3; %liquid only
            
            

        if length(igood2)>0
            
            
            if i_make_all ==1

                %all of the Nd data
                Nd_all(ilatL3L2,ilonL3L2,1:length(igood2)) = Nd_L2_swath(igood2);
                Tau_all(ilatL3L2,ilonL3L2,1:length(igood2)) = Tau_L2_swath(igood2);
                Re_all(ilatL3L2,ilonL3L2,1:length(igood2)) = Re_L2_swath(igood2);
                Nd_un_all(ilatL3L2,ilonL3L2,1:length(igood2)) = Nd_un_L2_swath(igood2);
                Re_un_all(ilatL3L2,ilonL3L2,1:length(igood2)) = Re_un_L2_swath(igood2);
                Tau_un_all(ilatL3L2,ilonL3L2,1:length(igood2)) = Tau_un_L2_swath(igood2);
                CTT_all(ilatL3L2,ilonL3L2,1:length(igood2)) = CTT_L2_swath(igood2);
                SolarZA_all(ilatL3L2,ilonL3L2,1:length(igood2)) = solarZA_L2_swath(igood2);                
                SensorZA_all(ilatL3L2,ilonL3L2,1:length(igood2)) = sensorZA_L2_swath(igood2);
            end
            
            
           
            
            
            Nreject_mockL3(ilatL3L2,ilonL3L2,itime) = length(ipos_L2)-length(igood2);
            

            

            %1D histograms within the cell
%             NdPDF_mockL3(ilatL3L2,ilonL3L2,:,itime) = ndhistc_run(Nd_L2_swath(igood2), Nd_bins);
%             NdPDF_mockL3_16(ilatL3L2,ilonL3L2,:,itime) = ndhistc_run(Nd_L2_swath_16(igood2), Nd_bins);
%             NdPDF_mockL3_37(ilatL3L2,ilonL3L2,:,itime) = ndhistc_run(Nd_L2_swath_37(igood2), Nd_bins);            
%             WPDF_mockL3(ilatL3L2,ilonL3L2,:,itime) = ndhistc_run(W_L2_swath(igood2), W_bins);  
%             WPDF_mockL3_16(ilatL3L2,ilonL3L2,:,itime) = ndhistc_run(W_L2_swath_16(igood2), W_bins);              
%             WPDF_mockL3_37(ilatL3L2,ilonL3L2,:,itime) = ndhistc_run(W_L2_swath_37(igood2), W_bins);              
%             RePDF_mockL3(ilatL3L2,ilonL3L2,:,itime) = ndhistc_run(Re_L2_swath(igood2), Re_bins);
%             RePDF_mockL3_16(ilatL3L2,ilonL3L2,:,itime) = ndhistc_run(Re_L2_swath_16(igood2), Re_bins);
%             RePDF_mockL3_37(ilatL3L2,ilonL3L2,:,itime) = ndhistc_run(Re_L2_swath_37(igood2), Re_bins);            
%             TauPDF_mockL3(ilatL3L2,ilonL3L2,:,itime) = ndhistc_run(Tau_L2_swath(igood2), Tau_bins);
%             
%             Tau_un_PDF_mockL3(ilatL3L2,ilonL3L2,:,itime) = ndhistc_run(Tau_un_L2_swath(igood2), Tau_un_bins);
%             Re_un_PDF_mockL3(ilatL3L2,ilonL3L2,:,itime) = ndhistc_run(Re_un_L2_swath(igood2), Re_un_bins);
%             Nd_un_PDF_mockL3(ilatL3L2,ilonL3L2,:,itime) = ndhistc_run(Nd_un_L2_swath(igood2), Nd_un_bins);

            % *** calculate a 1D histogram of Nd vs CTP.
            %first create a 2D histo of Nd and CTP
            %Nd_CTP_2D_PDF = ndhistc_run([Nd_L2_swath(igood2) CTP_L2_swath(igood2)], Nd_2D_bins, CTP_2D_bins);
            Nd_bins_rep2D = ( repmat(0.5*(Nd_2D_bins(1:end-1)+Nd_2D_bins(2:end)),[length(CTP_2D_bins)-1 1]) )';
            %Nd_vs_CTP_mockL3(ilatL3L2,ilonL3L2,:,itime) = sum(Nd_CTP_2D_PDF.*Nd_bins_rep2D,1)./sum(Nd_CTP_2D_PDF,1);
            %NP_vs_CTP_mockL3(ilatL3L2,ilonL3L2,:,itime) = sum(Nd_CTP_2D_PDF,1);            
            
            %2D tau-reff histograms
            %Tau_Re_2D_PDF(ilatL3L2,ilonL3L2,:,:,itime) = ndhistc_run([Tau_L2_swath(igood2) Re_L2_swath(igood2)], Tau2D_bins, Re2D_bins);
            


            [meanNd_mockL3(ilatL3L2,ilonL3L2,itime),Npoints_temp,stdNd_mockL3(ilatL3L2,ilonL3L2,itime)] = meanNoNan(Nd_L2_swath(igood2),1);
            [meanNd_mockL3_16(ilatL3L2,ilonL3L2,itime),Npoints_temp,stdNd_mockL3_16(ilatL3L2,ilonL3L2,itime)] = meanNoNan(Nd_L2_swath_16(igood2),1);
            [meanNd_mockL3_37(ilatL3L2,ilonL3L2,itime),Npoints_temp,stdNd_mockL3_37(ilatL3L2,ilonL3L2,itime)] = meanNoNan(Nd_L2_swath_37(igood2),1);            
            [meanTau_mockL3(ilatL3L2,ilonL3L2,itime),Npoints_temp,stdTau_mockL3(ilatL3L2,ilonL3L2,itime)]  = meanNoNan(Tau_L2_swath(igood2),1);
            [meanRe_mockL3(ilatL3L2,ilonL3L2,itime),Npoints_temp,stdRe_mockL3(ilatL3L2,ilonL3L2,itime)]  = meanNoNan(Re_L2_swath(igood2),1);
            [meanRe_mockL3_16(ilatL3L2,ilonL3L2,itime),Npoints_temp,stdRe_mockL3_16(ilatL3L2,ilonL3L2,itime)]  = meanNoNan(Re_L2_swath_16(igood2),1);
            [meanRe_mockL3_37(ilatL3L2,ilonL3L2,itime),Npoints_temp,stdRe_mockL3_37(ilatL3L2,ilonL3L2,itime)]  = meanNoNan(Re_L2_swath_37(igood2),1);                     
            minRe_mockL3(ilatL3L2,ilonL3L2,itime) = min(Re_L2_swath(igood2));
            minRe_mockL3_16(ilatL3L2,ilonL3L2,itime) = min(Re_L2_swath_16(igood2));
            minRe_mockL3_37(ilatL3L2,ilonL3L2,itime) = min(Re_L2_swath_37(igood2));    
            maxRe_mockL3(ilatL3L2,ilonL3L2,itime) = max(Re_L2_swath(igood2));
            maxRe_mockL3_16(ilatL3L2,ilonL3L2,itime) = max(Re_L2_swath_16(igood2));
            maxRe_mockL3_37(ilatL3L2,ilonL3L2,itime) = max(Re_L2_swath_37(igood2));               
            meanNd_un_mockL3(ilatL3L2,ilonL3L2,itime) = meanNoNan(Nd_un_L2_swath(igood2),1);
            
            %sum the squares of the absolute uncertainty in Nd - use
            %meanNoNan and then mutliply by N to get the non-Nan sum. Then
            %square root - using the propagation of error formlua dY =
            %sqrt( (dY/dX1*error_X1)^2 + ..... ). Here Y is the mean of N,
            %so Y = sum(Ni)/Nvals
            [Nmean,Nvals] = meanNoNan( ( Nd_un_L2_swath(igood2).*Nd_L2_swath(igood2)/100 ) .^2,1); 
            %= 1/Nvals .* sqrt(Nmean.*Nvals); Nmean.*Nvals is the sum and
            %*1/Nvals as we want the error in the mean N, so this results
            %in sqrt(Nmean/Nvals)
            meanNd_un_combined_mockL3(ilatL3L2,ilonL3L2,itime) = sqrt(Nmean./Nvals);
            %note - we don't have L2 uncertainties for 1.6 and 3.7 um
            
            meanTau_un_mockL3(ilatL3L2,ilonL3L2,itime) = meanNoNan(Tau_un_L2_swath(igood2),1);
            meanTau_logun_mockL3(ilatL3L2,ilonL3L2,itime) = meanNoNan(log(Tau_un_L2_swath(igood2)),1);            
            maxTau_mockL3(ilatL3L2,ilonL3L2,itime) = max(Tau_L2_swath(igood2));
            minTau_mockL3(ilatL3L2,ilonL3L2,itime) = min(Tau_L2_swath(igood2));            
            meanRe_un_mockL3(ilatL3L2,ilonL3L2,itime) = meanNoNan(Re_un_L2_swath(igood2),1);
            meanRe_logun_mockL3(ilatL3L2,ilonL3L2,itime) = meanNoNan(log(Re_un_L2_swath(igood2)),1);            
            [meanCTT_mockL3(ilatL3L2,ilonL3L2,itime),Npoints_temp,stdCTT_mockL3(ilatL3L2,ilonL3L2,itime)] = meanNoNan(CTT_L2_swath(igood2),1);
            minCTT_mockL3(ilatL3L2,ilonL3L2,itime) = min(CTT_L2_swath(igood2));
            maxCTT_mockL3(ilatL3L2,ilonL3L2,itime) = max(CTT_L2_swath(igood2));   
            
  
            
            
            
            
            [meanW_mockL3(ilatL3L2,ilonL3L2,itime),Npoints_temp,stdW_mockL3(ilatL3L2,ilonL3L2,itime)] = meanNoNan(W_L2_swath(igood2),1);
            %for Cahalan's homogeneity factor (= exp(mean(lnW)) ./ mean(W) )            
            meanlogW_mockL3(ilatL3L2,ilonL3L2,itime) = meanNoNan(log(W_L2_swath(igood2)),1);
            meanTau_log_mockL3(ilatL3L2,ilonL3L2,itime) = meanNoNan(log(Tau_L2_swath(igood2)),1);            
            meanRe_log_mockL3(ilatL3L2,ilonL3L2,itime) = meanNoNan(log(Re_L2_swath(igood2)),1);            
            meanRe_log_mockL3_16(ilatL3L2,ilonL3L2,itime) = meanNoNan(log(Re_L2_swath_16(igood2)),1);
            meanRe_log_mockL3_37(ilatL3L2,ilonL3L2,itime) = meanNoNan(log(Re_L2_swath_37(igood2)),1);
            
            %don't want to only use liquid phase etc for solar and sensor
            %data
%             meanSolarZA_mockL3(ilatL3L2,ilonL3L2,itime) = meanNoNan(solarZA_L2_swath(igood2),1);
%             meanSolarAz_mockL3(ilatL3L2,ilonL3L2,itime) = meanNoNan(solarAz_L2_swath(igood2),1);            
%             meanSensorZA_mockL3(ilatL3L2,ilonL3L2,itime) = meanNoNan(sensorZA_L2_swath(igood2),1);    
%             meanSensorAz_mockL3(ilatL3L2,ilonL3L2,itime) = meanNoNan(sensorAz_L2_swath(igood2),1);
%             

            
        end
        
        %these are variables that aren't restricted to being in liquid
        %cloud
        meanSolarZA_mockL3(ilatL3L2,ilonL3L2,itime) = meanNoNan(solarZA_L2_swath(ipos_L2),1);
        meanSolarAz_mockL3(ilatL3L2,ilonL3L2,itime) = meanNoNan(solarAz_L2_swath(ipos_L2),1);
        meanSensorZA_mockL3(ilatL3L2,ilonL3L2,itime) = meanNoNan(sensorZA_L2_swath(ipos_L2),1);
        meanSensorAz_mockL3(ilatL3L2,ilonL3L2,itime) = meanNoNan(sensorAz_L2_swath(ipos_L2),1);
        
        if length(ice_new)>0 %otherwise min and max produce an error (meanNoNaN just returns NaN)
        meanCTT_ice_mockL3(ilatL3L2,ilonL3L2,itime) = meanNoNan(CTT_L2_swath(ipos_L2(ice_new)),1);
        minCTT_ice_mockL3(ilatL3L2,ilonL3L2,itime) = min(CTT_L2_swath(ipos_L2(ice_new)));
        maxCTT_ice_mockL3(ilatL3L2,ilonL3L2,itime) = max(CTT_L2_swath(ipos_L2(ice_new)));
        end

        
       end %ilatL3L2
    
       
       
    switch box_type
        case 'NxN pixel square'

            %want the edges of the averaging boxes for plotting purposes -
            %should be the first pixel - or could do an average along the box
            %edge?
            %add in the last edges
            %    LAT_mockL3_edge(ilatL3L2,ilonL3L2+1,itime) = mean(Plat_L2(ix,iy(end)+1));
            %    LON_mockL3_edge(ilatL3L2,ilonL3L2+1,itime) = mean(Plon_L2(ix,iy(end)+1));

            LAT_mockL3_edge(ilatL3L2,ilonL3L2+1,itime) = Plat_L2(ix(1),iy(end)+1);
            LON_mockL3_edge(ilatL3L2,ilonL3L2+1,itime) = Plon_L2(ix(1),iy(end)+1);

    end
    
    

    end  %ilonL3L2 

switch box_type
    case 'NxN pixel square'

        %add another row
        for ilonL3L2=i_LONS  %thses are already 1 less than the max  %1:size(LAT_mockL3_edge,2)-1

            xpos=(ilatL3L2-1+1)*Npix+1;
            ypos=(ilonL3L2-1)*Npix+1;
            ix=[xpos];
            iy=[ypos:ypos+Npix-1];
            %     [ix2,iy2]=meshgrid(ix,iy);
            %     %linear indices
            %     ipos_L2 = sub2ind(size(Plat_L2),ix2(:),iy2(:));



            %    LAT_mockL3_edge(ilatL3L2+1,ilonL3L2,itime) = mean(Plat_L2(ipos_L2));
            %    LON_mockL3_edge(ilatL3L2+1,ilonL3L2,itime) = mean(Plon_L2(ipos_L2));

            LAT_mockL3_edge(ilatL3L2+1,ilonL3L2,itime) = Plat_L2(ix,iy(1));
            LON_mockL3_edge(ilatL3L2+1,ilonL3L2,itime) = Plon_L2(ix,iy(1));

        end
        %the corner value
        LAT_mockL3_edge(ilatL3L2+1,ilonL3L2+1,itime) = Plat_L2(ix,iy(end)+1);
        LON_mockL3_edge(ilatL3L2+1,ilonL3L2+1,itime) = Plon_L2(ix,iy(end)+1);

end


Nd_vs_CTP_mockL3_noNaN = Nd_vs_CTP_mockL3;
Nd_vs_CTP_mockL3_noNaN(isnan(Nd_vs_CTP_mockL3))=0;
NP_vs_CTP_mockL3(isnan(Nd_vs_CTP_mockL3))=0;
totNP_vs_CTP = sum( sum( NP_vs_CTP_mockL3(:,:,:,itime) ,1),2);
meanNd_vs_CTP_mockL3 = sum( sum ( Nd_vs_CTP_mockL3_noNaN(:,:,:,itime) .* NP_vs_CTP_mockL3(:,:,:,itime) ,1) ,2) ./ totNP_vs_CTP;
midCTP_2D_bins = 0.5* ( CTP_2D_bins(1:end-1)+CTP_2D_bins(2:end) );


%  Wflag='calc';
%  [N_all,H_all,W_all,k,Q,cw]=MODIS_N_H_func(Tau_all,Re_all*1e-6,Wflag,NaN,CTT_all);
 
 [N_meanTauReff,H_meanTauReff,W_meanTauReff,k,Q,cw_meanTauReff]=MODIS_N_H_func(meanTau_mockL3,meanRe_mockL3*1e-6,Wflag,NaN,meanCTT_mockL3);
 [N_meanTauReff_16,H_meanTauReff_16,W_meanTauReff_16,k,Q,cw_meanTauReff]=MODIS_N_H_func(meanTau_mockL3,meanRe_mockL3_16*1e-6,Wflag,NaN,meanCTT_mockL3);
 [N_meanTauReff_37,H_meanTauReff_37,W_meanTauReff_37,k,Q,cw_meanTauReff]=MODIS_N_H_func(meanTau_mockL3,meanRe_mockL3_37*1e-6,Wflag,NaN,meanCTT_mockL3);

% %mid-points of the bins
% NdPDF_mid=0.5*(Nd_bins(2:end)+Nd_bins(1:end-1));
% Ngood=sum(NdPDF_mockL3,3);
% 
% Nmax=maxALL(Ngood);
% Npdf_all = NaN*ones([size(Np_mockL3) Nmax]);
% 
% for i=1:size(Np_mockL3,1)
%     for j=1:size(Np_mockL3,2)
%         N=Ngood(i,j);
%         if N~=0 & ~isnan(N)
%             Npdf_all(i,j,1:N)=expand_PDF(NdPDF_mid,squeeze(NdPDF_mockL3(i,j,:)) );
%         end
%     end
% end

    
clear override_mockL3_options2
catch mock_L3L2_error
    clear override_mockL3_options2
    rethrow(mock_L3L2_error);
end


