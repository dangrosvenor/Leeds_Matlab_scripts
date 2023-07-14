clear  modis_data_case
switch data_type
    case 'L3 processed data'  %uses this also for processing the timeseries3 data into monthly averages
        
        imod=1;
        %selects files for load_saved_modis_vars.m & modis_make_monthly_averages_multi_year.m
        
        % N.B. - the s='Y' and 'N' below don't do anything - they were just for a
        % note of which ones I had processed.
        
        if ~exist('override_loadsave') | override_loadsave==0
            L3_case = 'choose';
            %L3_case = 'Aqua and Terra, all years';
            %L3_case = 'Aqua only, 2006-2010'; %Was previously using this, but realized (14th March, 2016) that
            % the CALIPSO data is actually only
            % 2007-2010 for full years (2006 is
            % only partial).
            
            %    L3_case = 'Aqua only, 2007-2010'; %Need to match teh SST, though...
            %    L3_case = 'Aqua only, 2005-2012';
            
            L3_case = 'Aqua only, C6.1, 2016, 2017';
        end
        
        switch L3_case
            case 'override'
                %Chosen outside of this script
                modis_data_case = modis_data_case_choose;
            case 'choose'
                
                s='Y'; modis_data_case{imod} = {'y2000_TERRA'}; imod=imod+1; %
                s='Y'; modis_data_case{imod} = {'y2001_TERRA'}; imod=imod+1; %
                
                s='Y'; modis_data_case{imod} = {'y2002_AQUA'}; imod=imod+1; %
                s='Y'; modis_data_case{imod} = {'y2002_TERRA'}; imod=imod+1; %
                s='Y'; modis_data_case{imod} = {'y2003_AQUA'}; imod=imod+1; %
                s='Y'; modis_data_case{imod} = {'y2003_TERRA'}; imod=imod+1; %
                
                s='Y'; modis_data_case{imod} = {'y2004_AQUA'}; imod=imod+1; %with max SZA, CTT, scattering angle'};
                s='Y'; modis_data_case{imod} = {'y2004_TERRA'}; imod=imod+1; %with max SZA, CTT, scattering angle'};
                
                s='Y'; modis_data_case{imod} = {'y2005_TERRA'}; imod=imod+1; %Y2005 TERRA with max SZA, CTT, scattering angle'};
                s='Y'; modis_data_case{imod} = {'y2005_AQUA'}; imod=imod+1; %Y2005 AQUA with max SZA, CTT, scattering angle'};
                
                clear modis_data_case
                imod=1;
                
                %       s='Y'; modis_data_case{imod} = {'y2006_TERRA'}; imod=imod+1; %Y2006 TERRA with max SZA, CTT, scattering angle'};
                %       s='Y'; modis_data_case{imod} = {'y2006_AQUA_updated'}; imod=imod+1; %Y2006 AQUA (full year-partial 5.1) with max SZA, CTT, scattering angle'};
                %         s='Y'; modis_data_case{imod} = {'y2007_TERRA'}; imod=imod+1; %Y2007 TERRA with max SZA, CTT, scattering angle'};
                s='Y'; modis_data_case{imod} = {'y2002_AQUA'}; imod=imod+1; %Y2007
                %    TERRA with max SZA, CTT, scattering angle'};
                %        s='Y'; modis_data_case{imod} = {'y2008_TERRA'}; imod=imod+1; %Y2008 TERRA with max SZA, CTT, scattering angle'};
                %note - days 356 & 357 missing from LAADS website for TERRA 2008 -
                %system test / failure?? So, only have 364 days (out of 366).
                %        s='Y'; modis_data_case{imod} = {'y2008_AQUA'}; imod=imod+1; %
                
                %        s='Y'; modis_data_case{imod} = {'y2009_TERRA'}; imod=imod+1; %
                %        s='N'; modis_data_case{imod} = {'y2009_AQUA'}; imod=imod+1; %waiting for download
                %        s='Y'; modis_data_case{imod} = {'y2010_TERRA'}; imod=imod+1; %
                %        s='Y'; modis_data_case{imod} = {'y2010_AQUA'}; imod=imod+1; %
                
                %        s='Y'; modis_data_case{imod} = {'y2014_TERRA'}; imod=imod+1; %
                %        s='Y'; modis_data_case{imod} = {'y2014_AQUA'}; imod=imod+1; %
                
                %   clear modis_data_case
                %      imod=1;
                %   s='Y'; modis_data_case{imod} = {'y2008_AQUA'}; imod=imod+1; %
                %  s='Y'; modis_data_case{imod} = {'y2008_TERRA'}; imod=imod+1; %
                %    clear modis_data_case
                %    imod=1;
                %       s='Y'; modis_data_case{imod} = {'y2003_AQUA'}; imod=imod+1; %
                %      s='Y'; modis_data_case{imod} = {'y2003_TERRA'}; imod=imod+1; %
                
            case 'Aqua and Terra, all years'
                
                clear modis_data_case
                imod=1;
                
                s='Y'; modis_data_case{imod} = {'y2006_TERRA'}; imod=imod+1; %Y2006 TERRA with max SZA, CTT, scattering angle'};
                s='Y'; modis_data_case{imod} = {'y2006_AQUA_updated'}; imod=imod+1; %Y2006 AQUA (full year-partial 5.1) with max SZA, CTT, scattering angle'};
                s='Y'; modis_data_case{imod} = {'y2007_TERRA'}; imod=imod+1; %Y2007 TERRA with max SZA, CTT, scattering angle'};
                s='Y'; modis_data_case{imod} = {'y2007_AQUA'}; imod=imod+1; %Y2007
                %    TERRA with max SZA, CTT, scattering angle'};
                s='Y'; modis_data_case{imod} = {'y2008_TERRA'}; imod=imod+1; %Y2008 TERRA with max SZA, CTT, scattering angle'};
                %note - days 356 & 357 missing from LAADS website for TERRA 2008 -
                %system test / failure?? So, only have 364 days (out of 366).
                s='Y'; modis_data_case{imod} = {'y2008_AQUA'}; imod=imod+1; %
                
                s='Y'; modis_data_case{imod} = {'y2009_TERRA'}; imod=imod+1; %
                s='N'; modis_data_case{imod} = {'y2009_AQUA'}; imod=imod+1; %waiting for download
                s='Y'; modis_data_case{imod} = {'y2010_TERRA'}; imod=imod+1; %
                s='Y'; modis_data_case{imod} = {'y2010_AQUA'}; imod=imod+1; %
                
                
            case 'Aqua only, 2006-2010' %e.g. for AMSRE/CALIPSO comparison
                
                
                clear modis_data_case
                imod=1;
                
                %        s='Y'; modis_data_case{imod} = {'y2006_TERRA'}; imod=imod+1; %Y2006 TERRA with max SZA, CTT, scattering angle'};
                s='Y'; modis_data_case{imod} = {'y2006_AQUA_updated'}; imod=imod+1; %Y2006 AQUA (full year-partial 5.1) with max SZA, CTT, scattering angle'};
                %        s='Y'; modis_data_case{imod} = {'y2007_TERRA'}; imod=imod+1; %Y2007 TERRA with max SZA, CTT, scattering angle'};
                s='Y'; modis_data_case{imod} = {'y2007_AQUA'}; imod=imod+1; %Y2007
                %    TERRA with max SZA, CTT, scattering angle'};
                %        s='Y'; modis_data_case{imod} = {'y2008_TERRA'}; imod=imod+1; %Y2008 TERRA with max SZA, CTT, scattering angle'};
                %note - days 356 & 357 missing from LAADS website for TERRA 2008 -
                %system test / failure?? So, only have 364 days (out of 366).
                s='Y'; modis_data_case{imod} = {'y2008_AQUA'}; imod=imod+1; %
                
                %        s='Y'; modis_data_case{imod} = {'y2009_TERRA'}; imod=imod+1; %
                s='N'; modis_data_case{imod} = {'y2009_AQUA'}; imod=imod+1; %waiting for download
                %        s='Y'; modis_data_case{imod} = {'y2010_TERRA'}; imod=imod+1; %
                s='Y'; modis_data_case{imod} = {'y2010_AQUA'}; imod=imod+1; %
                
            case 'Aqua only, 2007-2010' %e.g. for AMSRE/CALIPSO comparison
                
                
                clear modis_data_case
                imod=1;
                
                %        s='Y'; modis_data_case{imod} = {'y2006_TERRA'}; imod=imod+1; %Y2006 TERRA with max SZA, CTT, scattering angle'};
                %        s='Y'; modis_data_case{imod} = {'y2006_AQUA_updated'}; imod=imod+1; %Y2006 AQUA (full year-partial 5.1) with max SZA, CTT, scattering angle'};
                %        s='Y'; modis_data_case{imod} = {'y2007_TERRA'}; imod=imod+1; %Y2007 TERRA with max SZA, CTT, scattering angle'};
                s='Y'; modis_data_case{imod} = {'y2007_AQUA'}; imod=imod+1; %Y2007
                %    TERRA with max SZA, CTT, scattering angle'};
                %        s='Y'; modis_data_case{imod} = {'y2008_TERRA'}; imod=imod+1; %Y2008 TERRA with max SZA, CTT, scattering angle'};
                %note - days 356 & 357 missing from LAADS website for TERRA 2008 -
                %system test / failure?? So, only have 364 days (out of 366).
                s='Y'; modis_data_case{imod} = {'y2008_AQUA'}; imod=imod+1; %
                
                %        s='Y'; modis_data_case{imod} = {'y2009_TERRA'}; imod=imod+1; %
                s='N'; modis_data_case{imod} = {'y2009_AQUA'}; imod=imod+1; %waiting for download
                %        s='Y'; modis_data_case{imod} = {'y2010_TERRA'}; imod=imod+1; %
                s='Y'; modis_data_case{imod} = {'y2010_AQUA'}; imod=imod+1; %
                
                
            case 'Aqua only, 2005-2012' %e.g. for POLDER work
                
                clear modis_data_case
                imod=1;
                
                s='Y'; modis_data_case{imod} = {'y2005_AQUA'}; imod=imod+1; %Y2005 AQUA with max SZA, CTT, scattering angle'};
                %        s='Y'; modis_data_case{imod} = {'y2006_TERRA'}; imod=imod+1; %Y2006 TERRA with max SZA, CTT, scattering angle'};
                s='Y'; modis_data_case{imod} = {'y2006_AQUA_updated'}; imod=imod+1; %Y2006 AQUA (full year-partial 5.1) with max SZA, CTT, scattering angle'};
                %        s='Y'; modis_data_case{imod} = {'y2007_TERRA'}; imod=imod+1; %Y2007 TERRA with max SZA, CTT, scattering angle'};
                s='Y'; modis_data_case{imod} = {'y2007_AQUA'}; imod=imod+1; %Y2007
                %    TERRA with max SZA, CTT, scattering angle'};
                %        s='Y'; modis_data_case{imod} = {'y2008_TERRA'}; imod=imod+1; %Y2008 TERRA with max SZA, CTT, scattering angle'};
                %note - days 356 & 357 missing from LAADS website for TERRA 2008 -
                %system test / failure?? So, only have 364 days (out of 366).
                s='Y'; modis_data_case{imod} = {'y2008_AQUA'}; imod=imod+1; %
                
                %        s='Y'; modis_data_case{imod} = {'y2009_TERRA'}; imod=imod+1; %
                s='N'; modis_data_case{imod} = {'y2009_AQUA'}; imod=imod+1; %waiting for download
                %        s='Y'; modis_data_case{imod} = {'y2010_TERRA'}; imod=imod+1; %
                s='Y'; modis_data_case{imod} = {'y2010_AQUA'}; imod=imod+1; %
                s='Y'; modis_data_case{imod} = {'y2011_AQUA'}; imod=imod+1; %
                s='Y'; modis_data_case{imod} = {'y2012_AQUA'}; imod=imod+1; %
                
            case 'Aqua only, C6.1, 2016, 2017'
                
                clear modis_data_case
                imod=1;
                
                %N.B. - this corresponds to the directory.
                %s='Y'; modis_data_case{imod} = {'y2016_AQUA_C61'}; imod=imod+1; %Currently only Dec 2016
                s='Y'; modis_data_case{imod} = {'y2017_AQUA_C61'}; imod=imod+1; %Full year
                
                
        end
        
        
        
        %Have re-processed 2008 aqua and terra to use the CTT in the histogram Nd
        %calculations - but not for other years as yet (is this necessary anyway,
        %as can recalculate?
        
    case 'Monthly averages etc from timeseries3' %think just uses this for loading monthly data once processed
        
        imod=1;
        %selects files for load_saved_modis_vars.m & modis_make_monthly_averages_multi_year.m
        single_combined = 'combined';
        single_combined = 'single';
        
        switch single_combined
            case 'single' %multiple years' Nd_PDF_multi etc loaded together
                %these are just for loading data
                %     s='Y'; modis_data_case{imod} = {'y2000_TERRA'}; imod=imod+1; %
                %     s='Y'; modis_data_case{imod} = {'y2001_TERRA'}; imod=imod+1; %
                %     s='Y'; modis_data_case{imod} = {'y2002_TERRA'}; imod=imod+1; %
                %     s='Y'; modis_data_case{imod} = {'y2003_TERRA'}; imod=imod+1; %
                %     s='Y'; modis_data_case{imod} = {'y2004_TERRA'}; imod=imod+1; %
                %     s='Y'; modis_data_case{imod} = {'y2005_TERRA'}; imod=imod+1; %
                %     s='Y'; modis_data_case{imod} = {'y2006_TERRA'}; imod=imod+1; %
                
                s='Y'; modis_data_case{imod} = {'y2007_TERRA'}; imod=imod+1; %Y2007
                s='Y'; modis_data_case{imod} = {'y2008_TERRA'}; imod=imod+1; %
                
                s='Y'; modis_data_case{imod} = {'y2009_TERRA'}; imod=imod+1; %
                s='Y'; modis_data_case{imod} = {'y2010_TERRA'}; imod=imod+1; %
                
            case 'combined'
                %these are for the combined data - single files where multiple
                %years have been added together to save memory
                clear modis_data_case
                imod=1;
                %all of the years combined into one PDF to save memory for median
                s='Y'; modis_data_case{imod} = {'y2001-2010_AQUA_TERRA'}; imod=imod+1; %
        end
        
    case 'Mock L3 timeseries3' %
        
        imod=1;
        %selects files for load_saved_modis_vars.m & modis_make_monthly_averages_multi_year.m
        
        %    s='Y'; modis_data_case{imod} = {'AP_Dec_2009_daytime_AQUA'}; imod=imod+1; %
        
        %    s='Y'; modis_data_case{imod} = {'Summer_2004_17days_all_lats_aqua'}; imod=imod+1; %
        s='Y'; modis_data_case{imod} = {'Summer_2007_Arctic_20W-60E_70-80N_terra'}; imod=imod+1; %
        s='Y'; modis_data_case{imod} = {'Summer_2007_Arctic_20W-60E_70-80N_aqua'}; imod=imod+1; %
        
        
        %     clear modis_data_case
        %     imod=1;
        
        
end