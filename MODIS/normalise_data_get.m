try

    if ~exist('ioverride_normalise_flags') | ioverride_normalise_flags==0

        calc_or_load = 'load';

        normalise_case = 'CF';
        normalise_case = 'LWP';

        data_type_for_normalise = 'GCM';
        data_type_for_normalise = 'satellite';

    end



switch normalise_case
    case 'CF'
         MODIS_varname2_plot = [MODIS_varname2_plot ' normalised by CF'];
        
        switch data_type_for_normalise

            case 'GCM'

                normalise_data = eval(['cf_isccp_low_' gcm_str]);
                normalise_data(normalise_data<0.01)=NaN;
                ylab='Precip rate (mm hr^{-1}), normalised by CF ';


            case 'satellite'
                normalise_data2 = cllcalipso_monthly(13:60,:,:)/100;

                %want to average the lon dimension from 2 to 4 degrees
                %so make an array half the size in lon
                normalise_data3 = NaN*ones([size(normalise_data2,1),size(normalise_data2,2),size(normalise_data2,3)/2]);

                for inorm=1:2:size(normalise_data2,3)
                    ind_norm = (inorm)/2+0.5;
                    normalise_data3(:,:,ind_norm) = 0.5*(normalise_data2(:,:,inorm)+normalise_data2(:,:,inorm+1));

                end


                normalise_data3(normalise_data3<0.01)=NaN;

                normalise_data = NaN*ones([size(normalise_data2,1)*2,size(normalise_data2,2),size(normalise_data2,3)/2]);

                for inorm=1:size(normalise_data3,1)
                    ind_norm = (inorm-1)*2+1;
                    normalise_data(ind_norm,:,:) = normalise_data3(inorm,:,:);
                    normalise_data(ind_norm+1,:,:) = normalise_data3(inorm,:,:);
                end


        end

    case 'LWP'
          MODIS_varname2_plot = [MODIS_varname2_plot ' normalised by LWP '];
          
         switch data_type_for_normalise

             case 'GCM'
                 normalise_data = eval(['1e3*gcm_lwp_minthreshCF_' gcm_str ';']);
                 ylab='Precip rate (mm hr^{-1}), normalised by LWP';
                 
             case 'satellite'
               
                 units_str_plot='mm m^3 kg^{-1} hr^{-1}';
                 
                 switch calc_or_load
                     case 'calc'
                         
                 screen_type='NP + CF';
                 thresh_CF=[0.8 1.00001];
                 thresh_NP=50;

                 %screening done here --->
                 modisL3_screening_timeseries3
                 % **************************************
                 
                 %so, are only looking at aclouds with a MODIS CF>0.8.
                 %Normalisation data will be NaN outside of such regions.
                 %Need to do the same for the model preicp
                 
                 dat_modis = 1e3*W_time3.*Cloud_Fraction_Liquid.timeseries3*(1/1.15);
                 dat_modis(ihtot)=NaN; %remove unwanted data

                 [Y,modismonth_timseries3,D] = datevec( datenum(modisyear_timeseries3,1,1) + daynum_timeseries3_MODIS - 1 );
                 months_unique = unique(modismonth_timseries3);
%                 years_unique = unique(Y);  

%% Want to create a MODIS LWP array that is the same size as the CloudSat
%% precip array of [96 90 90] 96 is from 4 years of monthly data
%% (2007-2010) with both an ascending and descending value. Probably best
%% to use only Aqua data here since that will coincide with CloudSat

years_unique=[2007:2010];
                 ipos=1;
                 LWP_monthly = NaN*ones([size(dat_modis,1) size(dat_modis,2) 2*length(months_unique)*length(years_unique)]);
                 for iyear_unique=1:length(years_unique)      
                     year_now = years_unique(iyear_unique);
                     for imon_unique=1:length(months_unique)                        
                         month_now = months_unique(imon_unique);
                         imon = find(modismonth_timseries3==month_now & Y==year_now);
                         %the following will return NaNs if imon=[]
                         %ascending (daytime)
                         LWP_monthly(:,:,ipos) = meanNoNan(dat_modis(:,:,imon),3);
                         %descending (nighttime) - make NaN as we don't
                         %have  nightime LWP
                         LWP_monthly(:,:,ipos+1) = NaN*ones([size(dat_modis,1) size(dat_modis,2)]);                       
                         ipos = ipos + 2;
                     end
                 end
                 
                 LWP_monthly = permute(LWP_monthly,[3 1 2]);
                 
                 LWP_monthly_coarse_temp = NaN*ones([size(LWP_monthly,1) size(LWP_monthly,2) size(LWP_monthly,3)/4]);
                 
%%now average to the same lat lon grid as for CloudSat
                for inorm=1:4:size(LWP_monthly,3)
                    ind_norm = (inorm)/4+0.75;
                    LWP_monthly_coarse_temp(:,:,ind_norm) = meanNoNan(LWP_monthly(:,:,inorm:inorm+3),3);                 
                end
                
                LWP_monthly_coarse = NaN*ones([size(LWP_monthly,1) size(LWP_monthly,2)/2 size(LWP_monthly,3)/4]);                
               
                for inorm=1:2:size(LWP_monthly_coarse_temp,2)
                    ind_norm = (inorm)/2+0.5;
                    LWP_monthly_coarse(:,ind_norm,:) = meanNoNan(LWP_monthly_coarse_temp(:,inorm:inorm+1,:),2);
                end
                
                     case 'load'
                         loaddir = ['/home/disk/eos1/d.grosvenor/matlab/work/MODIS/cpt/'];
                         load([loaddir 'LWP_monthly_MODIS_2007-2008_CF_0.8-1.0.mat']);
                         
                 end


         end

end



clear ioverride_normalise_flags
catch error_normalise
    clear ioverride_normalise_flags
    rethrow(error_normalise);
end

    

