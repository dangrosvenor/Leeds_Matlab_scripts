%Takes en equivalent potemp value, matches it to an environmental sounding and then
%calculates the LWC that would be produced during a dry then moist ascent from that
%point (using the P, T and qv from the environmental sounding)

%run wave_temperature_field_estimate.m first to get the sounding profiles for flight 102

f=1e6*28.97/18; 
set_column_numbers_for_flight_data %sets col_alt, etc.

%%%% ****  get the environmental (ascent & descent) sounding ****
get_flight102_envir_sounding
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%get LWC for CAS and hotwire
CAS_LWC_cut_off_sizes=[0 50]; %upper cut off for CAS LWC (microns)

 time_timeseries = CAS_time_all;

airspeed_constant=0;

if exist('iuse_precomputed_lwc')==1 & iuse_precomputed_lwc==1
    display('**** WARNING - Using pre-computed LWC ****');
    clear iuse_precomputed_lwc
else        
    [sample_volume_CAS,sample_volume_CIP,air_speed_1D,air_speed,CAS_total_number(icas_count)...
        ,CAS_total_number_cutoff ...                                    
        ,CIP_total_number(icas_count),LWC_dist_cas,LWC_dist_cip,CAS_mode_diameter...
        ,CAS_mean_diameter,LWC_dist_cas_cutoff,LWC_size_dist,bin_range,LWC_dist_cas_cutoff2,MVD,MVD_cut_off]...
       =cas_sample_volume_and_stats2...
        (dat_flt,time_timeseries,...
        CAS_bins,CAS_counts_all,CIP_time_all,CIP_bins,CIP_counts_all,air_speed_type,cut_off_size,TAS_all,CAS_LWC_cut_off_sizes,airspeed_constant); 
end

LWC_cas = LWC_dist_cas_cutoff;
LWC_hot=interp1(CIP_time_all,LWC_CAS_all',CAS_time_all); %interpolate to be on the same time axis

ii=find(LWC_hot>0);
ii=1:length(LWC_hot);

vap_factor=1.025;

if vap_factor>1
    fprintf(1,'\n*** WARNING applying a vap_factor of %f ***\n',vap_factor);
end
%aircraft obs for each LWC point during the flight
T=273.15+interp1(dat_flt(:,1)/1e3,dat_flt(:,col_temp),CAS_time_all)';
P=100*interp1(dat_flt(:,1)/1e3,dat_flt(:,col_press),CAS_time_all)';
qv_humi=vap_factor*interp1(dat_flt(:,1)/1e3,qv_flt_humi,CAS_time_all)';
qv_fp=vap_factor*interp1(dat_flt(:,1)/1e3,qv_flt_humi,CAS_time_all)';
qv_sat = SatVapPress(T,'goff','liq',P,1)/f;


qv_instrument='fp';
qv_instrument='humi';

iassume_sat=1;

switch qv_instrument
    case 'fp'
        if iassume_sat==0
            qv=qv_fp;
        else            
            qv=qv_sat;
        end
        qv_sound = vap_factor*qfp_sound;
        
        

    case 'humi'
        if iassume_sat==0
            qv=qv_humi;
        else
            qv=qv_sat;
        end
        qv_sound = vap_factor*qhumi_sound;
end

equiv_sound = equivalent_potemp(273.15+Tsound,100*Psound,qv_sound);
equiv_flight = equivalent_potemp(T,P,qv);
 
clear LWC_diff LWC_end
for i=1:length(ii) %loop through all required LWC points (aircraft obs)
    
    %find point in the sounding with the closest equiv potemp
    is=findheight_nearest(equiv_sound,equiv_flight(ii(i)));
    P0=100*Psound(is);
    T0=273.15+Tsound(is);
    Q0=qv_sound(is);
    
    %calculate the LWC that would be produced from a dry then moist ascent from the sounding to the pressure
    %of the point in question
    [P_rad,qv_rad,T_rad,LWC_rad] = forward_dry_then_moist_adiabatic_path(P0,T0,Q0,P(ii(i)));
    if length(LWC_rad>0)
        LWC_end(i)=LWC_rad(end);       
    else
        LWC_end(i)=0;
    end
    
    LWC_diff_hot(i)=LWC_end(i)-LWC_hot(ii(i));
    LWC_diff_cas(i)=LWC_end(i)-LWC_cas(ii(i));
end

tol=0.02; %absolute error tolerance
iok_hot=find(abs(LWC_diff_hot)<tol);
iok_cas=find(abs(LWC_diff_cas)<tol);

tol=20; %percentage error tolerance
iok_hot2=find(abs(100*LWC_diff_hot./LWC_end)<tol);
iok_cas2=find(abs(100*LWC_diff_cas./LWC_end)<tol);

iplot=0;
if iplot==1
    plot(LWC_end,ydat(1).y,'kx');
    plot(LWC_end(iok_hot2),ydat(1).y(iok_hot2),'kx'); %careful - make sure ydat is correct
    %Hotwire
    plot(xdat(1).x(iok_hot2),ydat(1).y(iok_hot2),'kx'); %careful - make sure ydat is correct    
    %CAS
    plot(xdat(1).x(iok_cas2),ydat(1).y(iok_cas2),'kx'); %careful - make sure ydat is correct        
end

isave=0;
if isave==1
    append='';
    filename=remove_character(filename,'\','');
    filename=remove_character(filename,':','_');
    filename=remove_character(filename,'>','.GT.');
    filename(60:end)=''; %shorten
    saveas(gcf,[filedir filename append '.fig'],'fig');
%    print(gcf,[filedir filename append '.tiff'],'-dtiff','-r600');
    print(gcf,[filedir filename append '.jpeg'],'-djpeg','-r600');
end

    
disp('Finished Estimate_LWC')