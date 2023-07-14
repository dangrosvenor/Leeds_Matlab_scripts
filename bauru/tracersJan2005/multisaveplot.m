N=71; %max length of filename that we allow
         
comp='lacieLap';
comp='uni';
%comp='database4';

comp_save=comp;

% direc(1).dir='c:/cygwin/home/user/runs/dmidamp_2/';
% direc(2).dir='c:/cygwin/home/userw/runs/damp_ccn960/';
% %direc(4).dir='c:/documents and settings/tom/my documents/dan/dmi1715_5ppmv/';
% direc(3).dir='c:/cygwin/home/user/runs/damp_inx10/';
% direc(4).dir='g:rmuns/500mres/';

if ~exist('ichoose_itype_multisaveplot')
itype='prof'; %(or timeseries)
%itype='timh';
itype='scatter';
%itype='streamline';
%itype='';
else
    clear ichoose_itype_multisaveplot
end

icont_extra=0;

% ****
%idirs=[4:1:10]; %set to e.g. [1 1] for 2 E5MM plots for profiles)
idirs=[1];
% ****print
isavemulti=1;

%set iformats to save copies in several formats - e.g. set iformats={'ps','emf'} for postscript
%and emf format. Postscript is scalable. emf doesn't seem to be for complicated plots for some reason
%(but is ok for line plots!) why? postscript not good for viewing in word though as it doesn't diplay them
iformats={'ps','emf','fig'};
iformats={'emf','png'};
%iformats={'fig'};
%iformats={'jpeg'};
iformats={'png'};


isamescale2=0; %set this to zero if are doing a series of plots
subplotting=0; %flag for whether to use subplots instead of new figures
iabc=1; %flag to label figures (a), (b), etc (along with dirstamp)

if length(idirs)==1
    isamescale2=0;
end

idirsemm=[1 2 3]; %set for tim heights

lememm=0; %flag for lem emm comparisons  (set idirs to e.g. [1 1] - length of 2)
if lememm==1
    iplotselect=1;
end

if isamescale2==1 & strmatch(itype,'timh')==1
    reps=2;
else
    reps=1;
end



if subplotting==1
	bigcbar=1; %flag for one big colorbar instead of one for each subplot
else
    bigcbar=0; %flag for one big colorbar instead of one for each subplot
end
%bigcbar=0; 




savedirname='my documents/temp/';
%savedirname='my documents/Leeds_MMOCCA/EMM/';
%savedirname='my documents/logbook/vapour_paper/pics/Mirvette/';
%savedirname='my documents/logbook/vapour_paper/pics/supersats/'
%savedirname='my documents/logbook/MMOCCA/pics_CCN_1e-3/';
%savedirname='my documents/Leeds_MMOCCA/pics_lem_2008/';
savedirname='Y:\WRF\30thNov_Min\maxtot_profiles/total_water/';
savedirname='/home/mbexddg5/work/Ant_Jan06/vertical_wind_slices/';
%savedirname='my documents/WRF/ant_jan06_sfRUC_v3/equiv_potemp_cross-sections/';
savedirname='my documents/WRF/ant_jan06_ecmwf_ncep_polar/met_em_plots\wind_level3_572m_met_em/';
%savedirname='my documents/WRF/ant_jan06_sfRUC_v3/soil_temp_d03/';
savedirname='my documents\WRF\ecmwf_ml_0.5\d02_surface_pressure\';
savedirname='my documents\WRF\ncepR2_3dom_nudging_seaice\d02_surface_pressure\';
%savedirname='my documents\WRF\ant_jan06_ncep_polarWRF\equiv_cross_sections_y=295km_EC\';
%savedirname='my documents\WRF\aircraft_profiles\';
savedirname='my documents\WRF\plots\2ndFeb09\';
%savedirname='my documents\WRF\ecmwf_blres_test\wind_level6_273m\';
%savedirname='my documents\WRF\ncepR2_3dom_nudging\LH+SH\';
%savedirname='My Documents\WRF\ncepR2_seaice\d01_wind_level_15\';
%savedirname='My Documents\WRF\ecmwf_ml_0.5_nudging\wind_level_11_d02\';
savedirname='My Documents\WRF\ecmwf_ml_0.5_nudging\wind_10m_d03\';
%savedirname='My Documents\WRF\ecmwf_ml_0.5_nudging\wind_level_11_d03\';
%savedirname='My Documents\WRF\ncepR2_3dom_nudging_seaice\d01_level_15_wind\';
%savedirname='My Documents\WRF\Ant01_Dec_new_22Sep09\HM_zone_18UTC\';
%savedirname='My Documents\WRF\ncepR2_3dom_nudging_seaice\potemp_flight_cross_Rothera\';
savedirname='My Documents\WRF\ecmwf_ml_0.5_nudging\temperature_SKIN_d03\';
%savedirname='My Documents\WRF\ncepR2_3dom_nudging_seaice\temperature_level4\';
savedirname='My Documents\WRF\ecmwf_analysis_6th_Jan_case\temperature_06_d03_ANALYSIS_new\';
%savedirname='My Documents\WRF\ecmwf_analysis_6th_Jan_case\wind_all_levels_ecmwf_ANALYSIS_03_00_UTC\';



if ~exist('flight_no')
    flight_no='XXX';
end

comp_save='database4'; savedirname=['Y:\BAS_flights\flight' flight_no '\CAS_plots\'];
comp_save='database4'; savedirname=['Y:\BAS_flights\flight102\Wave_model_ACPIM_plots\'];

%comp_save='uni'; savedirname=['My Documents\logbook\Antarctica\Flights and instruments_Feb2010\Manchester and BAS CAS comparisons\'];
%comp_save='uni'; savedirname=['My Documents\logbook\Antarctica\Flights and instruments_Feb2010\Manchester and BAS CAS comparisons\29th_July_comparisons_with_FSSP_and_Welas\graphs\'];
%comp_save='database4'; savedirname='Y:\BAS_flights\calibrations\';
%comp_save='database4'; savedirname='Y:\BAS_flights\APPRAISE flight 433\';


%comp_save='uni'; savedirname='My Documents\WRF\ecmwf_ml_0.5_nudging\streamline_plots_6thJan_6UTC\';
%comp_save='uni'; savedirname='My Documents\WRF\ecmwf_ml_0.5_nudging\streamline_plots_6thJan_6UTC_black_coarse_2km_xy_closeup_terr_fine\';
%comp_save='uni'; savedirname='My Documents\WRF\ecmwf_ml_0.5_nudging\streamline_plots_6thJan_6UTC_2km_xy_closeup_14June2010_z_fine_dperp=400\';
%comp_save='uni'; savedirname='My Documents\WRF\ecmwf_ml_0.5_nudging\potemp_gradient_250-1500_d03\';


%comp_save='uni'; savedirname='My Documents\WRF\plots\29thSep09\';
%comp_save='uni'; savedirname='My Documents\WRF\plots\plots_for_AP_paper\';
%comp_save='uni'; savedirname='My Documents\logbook\Antarctica\Figures for paper Aug2010\';



clear h

tag='totbelow5_14-22km_14-Feb06';
%tag='_15.3-18km_29-Feb06';
tag='_13-22km_06-Mar06';
tag='_15.3-17km_06-Mar06';
tag='_0-22km_13-Mar06_24thFebcomp';
tag='28thMar06_vapour_250m_0cn_dump_1-31_13-22km';
tag='30thMarch_ccn_profiles';
tag='30thMarch06-250m5ppmv_1000km_dmi_5ppmv500km_comp_15-19km_dump1-49_noclines';
tag='pos_vels';
tag='_29-May06_250mres1000km_0-22km';
%tag='_NoNvapourPlume';
tag='250m_rhopert_1-250km';
%tag='MilesCity';

tag='15-30km_recomp_dump1-39';
tag='0-30km_recomp';
tag='CCNcomp_9.96km';
tag='all_ice_radcomps';
%tag=['_' num2str(it)];
tag='zoom';
tag='Texas_varying_strength2_dump14';
tag='streamline_JET_z0=800m';
%tag='pcolor_CORRECT_grid_x_grid_y';
%tag='LON=-66';
%tag='ART';
tag='LAT=67.5_z0=950m_129';
tag='';

minVal=9e99;
maxVal=-9e99;



iplotselect=0; %flag so that plotcase in plotTimeHeight is chosen from below
plotcases=[50 50 65 50];

iaxis_lims=1; %flag to say want to limit the axes by values below (automatically done for lememm==1)
timelims=[20.1 20.35];
zlims=[0 19]; %scale both plots to these values
clims=[0 1e8];
clims=[0 3.2];
clims=[0 14];
clims=[0 2.6];



xsub=2;
ysub=1;
onexlabel=1;


isamescale=0; %don't change this from zero
for irep=1:reps %repeat whole process if 
    clear h
    if isamescale2==1 & irep==2; isamescale=1; end

i577s={'vap_3d_vert' , 'vapour'};
i577s={'vap_3d_vert' , 'vap_3d_vert', 'tot_condensate'};
%i577s={'vap_3d_vert' , 'potemp'};
    
i3ds=[1 0 0 0];
i2ds=[1 1 1 1];

nsub=0;    
for iplot=1:length(idirs)
    
  %  i577=i577s{iplot};
  irad=iplot;
  
  try
      i3d=i3ds(iplot);
      wrap2d=i2ds(iplot);
  catch
  end
  
  
    
    
    idir=idirs(iplot);
    nsub=nsub+1;
    switch itype
        case 'timh'
            %wrap_slice;
            'plotting....'
            plottimeheightvap3;
            %  if isamescale2==1 & irep==1 & subplotting==1; close(gcf); end %close figure as will be drawing again
            if isamescale2==1 & irep==1 & (subplotting==0 | (subplotting==1 & iplot==length(idirs)) ); close(gcf); end %close figure as will be drawing again
        case 'prof'
            waterVapourMay2005;
        case 'scatter'
            scatter_plot;
        case 'streamline'
            iset_zstart=1;
            zstart=400 + 50*(idir-1); %increment in steps of 50m
            streamlines_threeD_draw
    end
    
         itext=findstr(':',savename);
         savename(itext)=',';
         
         itext=findstr('/',savename);
         savename(itext)='_';
         
         itext=findstr('\',savename);
         savename(itext)='_';
         
         itext=findstr('..',savename);
         savename(itext)='_';
         

         if length(savename)>N
            savename=savename(1:N);
          end

%'plotted?'
%pause(18);
%'pause over'
         

comp=comp_save;
    
    
     

     
 
 
    
   %  picname=[direcDan(idir).dir 'results/' savename '_' tag '.emf'];
     set(gcf,'paperpositionmode','auto');
     
     for icopy=1:length(iformats)
         iformat=iformats{icopy};
         
         
         % '&' characters are ok
         ichar=findstr(savename,'\');
         savename(ichar)='';
         ichar=findstr(savename,'>');
         savename(ichar)='';
         
         switch comp
             case 'uni'
                 picname=['c:/documents and settings/dan/' savedirname num2str(idir) '_' savename '_' tag '.' iformat];
             case 'lacieLap'
                 picname=['c:/documents and settings/g/' savedirname savename '_' num2str(idir) '_' tag '.' iformat];
             case 'database4'
                 picname=[savedirname num2str(idir) '_' savename '_' tag '.' iformat];

         end
    
         
         
         
     
     if isavemulti==1 & ( (isamescale==1 & subplotting==0) | (isamescale2==0 & subplotting==0) | ( isamescale==1 & subplotting==1 & iplot==length(idirs) ) ...
             | ( strmatch(itype,'prof')==1 & (subplotting==0 | (subplotting==1 & iplot==length(idirs)) ) )   ...
             | (subplotting==0 | (subplotting==1 & iplot==length(idirs) ) ) )
        %print(gcf,picname,'-dmeta');

		switch iformat
        case 'ps'
			print(gcf,picname,'-depsc');
        case 'jpeg'
			print(gcf,picname,'-djpeg75');    
        case 'png'
            print(gcf,picname,'-dpng','-r600'); %gives good looking output for all figures
            %including "patch" style figures
        case 'emf'
%			print(gcf,picname,'-dmeta');
            print(gcf,picname,'-dmeta');
            
            isave_fig=0;
            if isave_fig==1
                figname=['c:/documents and settings/login/'  savedirname num2str(iplot) '_' savename '_2' tag '_fig'];
                %remove brackets as doesn't allow filenames with brackets to be loaded
                ibrack=findstr(figname,'('); figname(ibrack)=' ';
                ibrack=findstr(figname,')'); figname(ibrack)=' ';

                saveas(gcf,figname,'fig');
            end
            
        
%            saveas(gcf,['c:/documents and settings/login/' savedirname savename '_' num2str(iplot) '_2' tag '_fig'],'mmat');        

            case 'fig'
                saveas(gcf,[picname],'fig');        
        end
     	
%pause(60); %pause for a while so can draw image - no need for pausing - problem was with isamescale2 not being set to zero so was closing figures before was saving them
        if icopy==length(iformats)
            close(gcf);
        end

     end
     
     end
   
   if strmatch(itype,'timh','exact')  
     if dlogflag==1
         minVal=min([minVal idlog(conts(1),dlogmin)]);
         maxVal=max([maxVal idlog(conts(end),dlogmin)]);
     else
            minVal=min([minVal conts(1)]);
            maxVal=max([maxVal conts(end)]);
     end    
   end
   
 end
 
 if strmatch(itype,'timh','exact')
     
     if dlogflag==1
         fprintf(1,'mincovOvr = dlog(%f,dlogmin);\nmaxcovOvr = dlog(%f,dlogmin);',minVal,maxVal);
     else
         fprintf(1,'mincovOvr = %f;\nmaxcovOvr = %f;',minVal,maxVal);
     end
 
end

end
         
isamescale=0; %reset this flag
