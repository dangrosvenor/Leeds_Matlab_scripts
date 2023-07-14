%load POLDER data that is colocated with MODIS for VOLCALS 2007
%N.B. POLDER is onboard the AQUA satellite



load_case ='VOCALS 2007'; %Note that the POLDER and MODIS data are averaged into daily values - so there is no individual orbit data.
    % But this should not be a problem for VOCALS where we are at less than 23
    % degrees - i.e. there should be no orbit overlap and so only one overpass
    % per day
load_case ='VOCALS 2005-2012';

switch load_case
     case 'VOCALS 2005-2012'
        %VOCALS region with the POLDER data co-located onto MODIS data
        %(mock L3)
        
        %No confidence screening 2005-2012
        loadfile = '/home/disk/eos8/d.grosvenor/mat_files_various/CPT/AMSRE_and_MODIS_2005-2012_Aqua_mockL3_no_confidence_screening_with_POLDER.mat';
%        load(loadfile);
        mockL3_no_conf = load(loadfile);
        
        %Load some extra files into different structures
        
        %No confidence screening 2007 only?
        %Case with no screening done using the water path confidence flags,
        %since it causes a truncation at re=20um (for re2.1)
%        loadfile_extra = '/home/disk/eos8/d.grosvenor/mat_files_various/CPT/VOCALS_MODIS_Aqua_2007_no_confidence_screening.mat';
%        mockL3_no_conf = load(loadfile_extra);
        


    case 'VOCALS 2007'
        %VOCALS region for 2007 with the POLDER data co-located onto MODIS data
        %(mock L3) - with confidence screening
        loadfile = '/home/disk/eos8/d.grosvenor/mat_files_various/CPT/AMSRE_and_MODIS_2006-2010_Aqua_with_reff37_etc_for2007_only.mat';
        load(loadfile)
        
        %Load some extra files into different structures
        
        %Case with no screening done using the water path confidence flags,
        %since it causes a truncation at re=20um (for re2.1)
        loadfile_extra = '/home/disk/eos8/d.grosvenor/mat_files_various/CPT/VOCALS_MODIS_Aqua_2007_no_confidence_screening.mat';
        mockL3_no_conf = load(loadfile_extra);
        
        %Case where the water path confidence screening has been done, but
        %the re<20um restriction has also been applied to re1.6 and re3.7
%        loadfile_extra = '/home/disk/eos8/d.grosvenor/mat_files_various/CPT/VOCALS_MODIS_Aqua_2007_confidence_screening_re16_re37_LT_20um_only.mat';
%        mockL3_conf_re16_re37_LT_20um = load(loadfile_extra);
           
end




%modisyear_timeseries3_MODIS = modisyear_timeseries3;

if exist('mockL3_no_conf')
    MODIS.LAT = mockL3_no_conf.MLAT;
    MODIS.LON = mockL3_no_conf.MLON;
    mockL3_no_conf.LAT = mockL3_no_conf.MLAT;
    mockL3_no_conf.LON = mockL3_no_conf.MLON;
    mockL3_no_conf.daynum_timeseries3 = mockL3_no_conf.daynum_timeseries3_MODIS;    
    mockL3_no_conf.modisyear_timeseries3 = mockL3_no_conf.modisyear_timeseries3_MODIS;    
    daynum_timeseries3_MODIS = mockL3_no_conf.daynum_timeseries3_MODIS;  
    modisyear_timeseries3_MODIS = mockL3_no_conf.modisyear_timeseries3_MODIS;
end

if exist('LAT') & exist('LAT')
    MODIS.LAT = LAT;
    MODIS.LON = LON;
elseif exist('MLAT')
    MODIS.LAT = MLAT;
    MODIS.LON = MLON;
end

LAT_MODIS = MODIS.LAT;
LON_MODIS = MODIS.LON;


% gcm_Plat2D_AMSRE = Plat2D_AMSRE_time3;
%     gcm_Plon2D_AMSRE = Plon2D_AMSRE_time3;    
%     gcm_Plat2D_edges_AMSRE = Plat2D_AMSRE_time3_edges;
%     gcm_Plon2D_edges_AMSRE = Plon2D_AMSRE_time3_edges;    
%     LAT_AMSRE = minALL(Plat2D_AMSRE_time3) :1: maxALL(Plat2D_AMSRE_time3);
%     LON_AMSRE = minALL(Plon2D_AMSRE_time3) :1: maxALL(Plon2D_AMSRE_time3);    
% 
% 
% if exist('daynum_timeseries3')
%     daynum_timeseries3_time3 = daynum_timeseries3;
% end
% 
% if ~exist('daynum_timeseries3_AMSRE') | ~ exist('gcm_time_UTC_AMSRE') 
%     nT = size(lwp_amsre,3);
%     daynum_timeseries3_AMSRE = 1:nT;
%     gcm_time_UTC_AMSRE=zeros([1 nT]);
% end


