      
        iadd_RWP=1; %whether to include model RWP as well as LWP
        iadd_Conv_LWP_RWP=1; %whether to include model LWP and RWP from convection scheme
        
        iplot_Ndat=0; %Whether to plot the number of days that pass the filter criteria, or the LWP
        ilow_clear_only=1; %Whether to restrict the MODEL to low+clear only scenes, or to all the CF>thresh choice to apply when there is also
        %mid or high cloud
        imodis_cf = 1; %Whether to use the COSP MODIS liquid cloud fractions (low, mid and high), or the model ones.
        iNtot_ratio_CF80_to_0 = 1; %flag that sets method for calculating the ration of Sc days to the total. Setting to one
        %uses the ratio of CF>80 to CF>=0 clouds. Clouds are still restricted
        %to low alt clouds by using mean CTH<=3.2km for MODIS and using only
        %gridboxees with no mid or high level cloud in the model (using COSP
        %MODIS low, mid and high clouds).. Should try to make them more
        %consistent.
        
        LWP_sat = 'AMSR-E';
        %LWP_sat = 'MODIS';
        
        %isave_plot_LWP_bias=0;
        fsize_latlon = 14;
 