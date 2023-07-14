%make monthly averages out of several years worth of data using the
%timeseries3 data as saved in different .mat files (one for each year,
%separately for aqua and terra)

%Prodcedure will be to 
%(1)   Loop through all of the desired years (for aqua and
%terra) loading in the timseries3 variables for each one at a time using
%load_saved_modis_vars.m 
%(2)   Then calculate monthly means of Nd, H and LWP and
%also Ndays (no. distinct days used) and Ndata (no. datapoints - counting
%the same day twice if viewed by Aqua and Terra or on multiple overpasses).
%Also will need to store Nd.^2 to be able to calculate std devs over
%multiple years.
%Also will do various degrees of screening for this - start with one file for CF>0.8 (can do 
%others later in separate files). Produce a matrix for the different
%degrees of screening.
%CTP-std_dev >=880,830,780,730,680 then also CTP+std_dev <= 680 & CTP-std_dev > 440  (ISCPP
%classifications) and                  CTP+std_dev <= 440 & CTP-std_dev > 50.
%Do for std_dev=0 (i.e. just based on mean CTP) and std_dev and 2*std_dev
%Gives 7*3=21 bins for pressure.
%12 bins for month
%could also do for sensorZA - say 7 bins [0:15:60]
%gives 21*12*5 = 1260 elements in total - too many? Could remove sensorZA
%screening? Or save in separate files?



 data_type='L3 processed data';
 mod_data_type='timeseries3 lambert';


% *** external script to select the required files :-
    savedir_var='/home/disk/eos8/d.grosvenor/saved_data_L3/';
    modisyear_timeseries3 = [];
    daynum_timeseries3 = [];
    aqua_terra_timeseries3 = [];
    
    modis_select_timeseries3_files
% **************************************************** 



thresh_CF=[0.8 1.0001];
thresh_NP=50;



%thresh_CTPs_1 = [880];
%thresh_CTPs_2 = [1000];



%bins for 1D PDFs - one for each variable in mdp_str - name as
%mdp_str_bins_mon
Nd_bins_mon = [0:25:200 300:100:1400 1500:500:5000];  %29 bins
LWP_bins_mon = [0:25:500 500:500:2000];  %29 bins
H_bins_mon = [0:25:500 500:500:2000];  %29 bins


thresh_ndays=0; %threshold no. days
proj_type='global oval';
data_select='specific_modis_data';
ifilter_ndays=1; %flag for whether to cut out points for which there aren't many days
thresh_ndays=0;
thresh_Ngood = 0;
icolormap_cf_grey=0;
irestrict_domain=0;


%  this funciton chooses thresh_CTP and the screen types
modis_CTP_screens




            


%parse the filenames to group the years together (will load aqua and terra
%in together)
modis_year_old=9999;
iyear=0;
nyear=1;
clear modis_data_case3
for imod=1:length(modis_data_case)  %loop through the files to load from
    filename_choose_saved_MODIS_files %basically using to determine the year
           
    if modis_year==modis_year_old | imod==1        
        iyear=iyear+1;
    else  %reset iyear & increment nyear
        nyear=nyear+1;        
        iyear=1;
    end
    
    modis_data_case3{nyear,iyear} = modis_data_case{imod};
    
    modis_year_old=modis_year;
    
    
end

%loop over all of the files selected
for iy_multi=1:nyear
    %clear/create arrays - will save each year separately
    
    %give the names of the variables to store
    modis_make_monthly_avs_saveload_variables
    %create the arrays
    for imdp2=1:imdp-1
        eval_str = [mdp_str{imdp2} '_multi = zeros([180 360 12 length(thresh_CTPs_1) length(screen_types)]);']; eval(eval_str);
        eval_str = [mdp_str{imdp2} '_sq_multi = zeros([180 360 12 length(thresh_CTPs_1) length(screen_types)]);']; eval(eval_str);
        eval_str = [mdp_str{imdp2} '_PDF_multi = zeros([180 360 12 length(thresh_CTPs_1) length(screen_types) length(' mdp_str{imdp2} '_bins_mon)-1]);']; eval(eval_str);
    end

    eval_str = ['Ndays_multi = zeros([180 360 12 length(thresh_CTPs_1) length(screen_types)]);']; eval(eval_str);
    eval_str = ['Ndata_multi = zeros([180 360 12 length(thresh_CTPs_1) length(screen_types)]);']; eval(eval_str);


    
%% load in the timeseries3 data - one array element for each day of the year

    %select the file to load based on modis_data_case3
    modis_case_savename = '';
    for im3=1:size(modis_data_case3,2)
        modis_data_case2{im3}=modis_data_case3{iy_multi,im3};
        modis_case_savename=[modis_case_savename '_' char(modis_data_case3{iy_multi,im3})];
    end

    i_multi_year_av=1; %flag for load_saved_modis_vars.m
    data_type='L3 processed data';
    override_loadsave=1;
    
% ----------------- Load modis data ---------------------------------
    load_saved_modis_vars
    %loads the case that is now contained in modis_data_case2 - N.B. if
    %both aqua and terra have been selected then they will both be
    %contained in modis_data_case2 and both loaded together
% -------------------------------------------------------------------
    


    
     for month_no=1:12
         month_no

            switch month_no
                case 1
                    days_required_for_mean = [1:31]; time_mean_str = 'Jan';
                case 2
                    days_required_for_mean = [32:60]; time_mean_str = 'Feb';
                case 3
                    days_required_for_mean = [61:91]; time_mean_str = 'Mar';
                case 4
                    days_required_for_mean = [92:121]; time_mean_str = 'Apr';
                case 5
                    days_required_for_mean = [122:152]; time_mean_str = 'May';
                case 6
                    days_required_for_mean = [153:182]; time_mean_str = 'Jun';
                case 7
                    days_required_for_mean = [183:213]; time_mean_str = 'Jul';
                case 8
                    days_required_for_mean = [214:244]; time_mean_str = 'Aug';
                case 9
                    days_required_for_mean = [245:274]; time_mean_str = 'Sep';
                case 10
                    days_required_for_mean = [275:305]; time_mean_str = 'Oct';
                case 11
                    days_required_for_mean = [306:335]; time_mean_str = 'Nov';
                case 12
                    days_required_for_mean = [336:366]; time_mean_str = 'Dec';
            end
            

            

            
            for iCTP_screens = 1:length(screen_types)
                screen_type=screen_types{iCTP_screens};                

                for iCTP=1:length(thresh_CTPs_1)
                    thresh_CTP(1)=thresh_CTPs_1(iCTP);
                    thresh_CTP(2)=thresh_CTPs_2(iCTP);


                    
                    for imdp=1:length(modis_data_plot_strs)
                        ioverride_plotglobal_thresh=1;
                        ioverride_time_selection=1;
                        noplot=1;
                        
                        modis_data_plot = modis_data_plot_strs{imdp};
                        plot_global_maps
    
                        %add the number of (unique) days to the running total
                        Ndays_multi(:,:,month_no,iCTP,iCTP_screens) = Ndays_multi(:,:,month_no,iCTP,iCTP_screens) + Ndays2;
                        %now for the number of actual data points
                        Ndata = ones(size(all_days));
                        Ndata(ihtot)=0;
                        Ndata_multi(:,:,month_no,iCTP,iCTP_screens) = Ndata_multi(:,:,month_no,iCTP,iCTP_screens) + sum(Ndata(:,:,time_inds_average),3);
                        
                        P2=P_save;
                        P2(isnan(P_save))=0; %remove NaNs - NaNs will not count to Ndata, so we can zero them
                      
%                        for imdp=1:length(imdp)
                            eval_str = [mdp_str{imdp} '_multi(:,:,month_no,iCTP,iCTP_screens) = ' mdp_str{imdp} '_multi(:,:,month_no,iCTP,iCTP_screens) + P2;']; eval(eval_str);
                            %the mean of the sqaures (individual Nd squared
                            %then mean of those) - different to P.^2 =
                            % ( mean(x_i=1:N) ).^2
                            eval_str = [mdp_str{imdp} '_sq_multi(:,:,month_no,iCTP,iCTP_screens) = ' mdp_str{imdp} '_sq_multi(:,:,month_no,iCTP,iCTP_screens) + meanNoNan(dat_modis(:,:,time_inds_average).^2,3);']; eval(eval_str);

                            %now make a PDF of the variable
                            %bins required (stated above for each variable)
                            eval_str = ['Xbins=' mdp_str{imdp} '_bins_mon;']; eval(eval_str);
                            
                            for ilat_mon=1:size(dat_modis,1)  %consider coding arrays into the C code for ndhistc?
                                for ilon_mon=1:size(dat_modis,2)


                                    X=dat_modis(ilat_mon,ilon_mon,time_inds_average);
                                    X(isnan(X))=[];
                                    X=X(:);
                                    qh = ndHistc_run([X], Xbins);
                                    
                                    eval_str = [mdp_str{imdp} '_PDF_multi(ilat_mon,ilon_mon,month_no,iCTP,iCTP_screens,1:length(Xbins)-1) = qh;']; eval(eval_str);



                                end
                            end


                            
                            
                            
%                        end
                        
                        close(gcf);
                        
                       
                        
                        
      
                       

                    end

                    


                end
            end
            
            
     end
     
    

     save_filename = [savedir_var 'monthly_mean_timeseries3_data_' modis_case_savename '_' datestr(now,30)];
     
     %script
     modis_make_monthly_avs_saveload_variables  %script
    
     fprintf(1,'\nSaving......');
     %save modis_var{1}
     eval_str = ['save(save_filename,''' modis_var{1} ''',''-v7.3'');']; eval(eval_str);
     %save the others, but appending this time
     for isave=2:length(modis_var)
         eval_str = ['save(save_filename,''' modis_var{isave} ''',''-APPEND'',''-v7.3'');']; eval(eval_str);
     end
     fprintf(1,'done\n');     
     
     
     
    
    
end %end of year loop