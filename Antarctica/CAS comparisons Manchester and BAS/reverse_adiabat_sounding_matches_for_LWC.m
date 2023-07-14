%assesing the consistency of a given observed LWC, P and T
%with an environmental (out of cloud) sounding.


isave=0;


%get LWC for CAS and hotwire
CAS_LWC_cut_off_sizes=[0 50]; %upper cut off for CAS LWC (microns)

 time_timeseries = CAS_time_all;
            
[sample_volume_CAS,sample_volume_CIP,air_speed_1D,air_speed,CAS_total_number(icas_count)...
    ,CAS_total_number_cutoff ...                                    
    ,CIP_total_number(icas_count),LWC_dist_cas,LWC_dist_cip,CAS_mode_diameter...
    ,CAS_mean_diameter,LWC_dist_cas_cutoff,LWC_size_dist,bin_range,LWC_dist_cas_cutoff2,MVD,MVD_cut_off]...
   =cas_sample_volume_and_stats2...
    (dat_flt,time_timeseries,...
    CAS_bins,CAS_counts_all,CIP_time_all,CIP_bins,CIP_counts_all,air_speed_type,cut_off_size,TAS_all,CAS_LWC_cut_off_sizes,airspeed_constant); 

LWC_cas = LWC_dist_cas_cutoff;
LWC_hot=interp1(CIP_time_all,LWC_CAS_all',CAS_time_all); %interpolate to be on the same time axis


%get aircraft temperatures, etc.
height_LWC = interp1(dat_flt(:,1)/1e3,dat_flt(:,col_alt),CAS_time_all)';
T=273.15+interp1(dat_flt(:,1)/1e3,dat_flt(:,col_temp),CAS_time_all)'; 
P=100*interp1(dat_flt(:,1)/1e3,dat_flt(:,col_press),CAS_time_all)';
qv = SatVapPress(T,'goff','liq',P,1)/f;   
   

                        
%do hotwire points
ii=find(LWC_hot>0.27 & LWC_hot<0.32); %find LWCs within this range
ii2=find(LWC_cas>0.54 & LWC_cas<0.6); %find LWCs within this range
ii2=ii; %make the same as for the hotwire points
%ii=ii2; %make the same as for the cas points

ii=find(LWC_hot>0.22 & LWC_hot<0.25 & ydat(1).y'>3.33);
ii=find(LWC_hot>0.18 & LWC_hot<0.23 & ydat(1).y'>3.3 & ydat(1).y'<3.34);
ii2=ii; %make the same as for the hotwire points

for i=1:length(ii) %loop through all the points selected and find the pressure, temperature and qv along 
                   %moist then dry adiabatic path (hotwire points)
    [P_rad_hot(i).dat,qv_rad_hot(i).dat,T_rad_hot(i).dat,LWC_rad_hot(i).dat,LWC_rad_kg_hot(i).dat] = ...
        reverse_moist_adiabatic_path(P(ii(i)),T(ii(i)),qv(ii(i)),LWC_hot(ii(i)));        
end



for i=1:length(ii2) %same for CAS points
    [P_rad_cas(i).dat,qv_rad_cas(i).dat,T_rad_cas(i).dat,LWC_rad_cas(i).dat,LWC_rad_kg_cas(i).dat] = ...
        reverse_moist_adiabatic_path(P(ii2(i)),T(ii2(i)),qv(ii2(i)),LWC_cas(ii2(i)));        
end


%plot for the vapour
figure
plot(qhumi_sound,Psound,'r'); hold on
plot(qfp_sound,Psound,'r--'); hold on
grid
set(gca,'ylim',[650 750]);
set(gca,'ydir','reverse');

for i=1:length(ii)
    plot(qv_rad_hot(i).dat,P_rad_hot(i).dat/100,'b');
end

for i=1:length(ii2)
    plot(qv_rad_cas(i).dat,P_rad_cas(i).dat/100,'g');
end

xlab='Vapour mixing ratio (kg kg^{-1})';
xlabel(xlab);
append='';
ylabel('Pressure (hPa)');

if isave==1
    saveas(gcf,[filedir xlab append '.fig'],'fig');
    print(gcf,[filedir xlab append '.emf'],'-dmeta'); 
    print(gcf,[filedir xlab append '.eps'],'-depsc');   
    fix_lines([filedir xlab append '.eps'],[filedir xlab append '.eps']);
end

%plot temperature
figure
plot(Tsound,Psound,'r'); hold on
grid
set(gca,'ylim',[650 750]);
set(gca,'ydir','reverse');

for i=1:length(ii)
    plot(T_rad_hot(i).dat,P_rad_hot(i).dat/100,'b');
end

for i=1:length(ii2)
    plot(T_rad_cas(i).dat,P_rad_cas(i).dat/100,'g');
end

xlab='Temperature (^{o}C)';
xlabel(xlab);
ylabel('Pressure (hPa)');

if isave==1
    saveas(gcf,[filedir xlab append '.fig'],'fig');
    print(gcf,[filedir xlab append '.emf'],'-dmeta'); 
    print(gcf,[filedir xlab append '.eps'],'-depsc');   
    fix_lines([filedir xlab append '.eps'],[filedir xlab append '.eps']);
end

scatter_plot
plot(LWC_hot(ii),ydat(1).y(ii),'kx');
xlab='Scatter plot with points marked';
if isave==1
    saveas(gcf,[filedir xlab append '.fig'],'fig');
    print(gcf,[filedir xlab append '.emf'],'-dmeta'); 
    print(gcf,[filedir xlab append '.eps'],'-depsc');   
    fix_lines([filedir xlab append '.eps'],[filedir xlab append '.eps']);
end















