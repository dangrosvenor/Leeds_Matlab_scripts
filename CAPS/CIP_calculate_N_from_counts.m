%calculate the total number concentrations from the CIP counts
%using the size dependent sample volume

tas_type='actual';
%tas_type='constant';

switch tas_type
    case 'actual'
        TAS_sVol = TAS_all; %use the actual airspeed
    case 'constant'
    TAS_sVol = 60*ones(size(TAS_all)); %use constant airspeed to be consistent with Jonny
end

%this is done in cas_sample_volume_and_stats2.m already (=sample_volume_CIP)
clear sVol_CIP
for i=0:60  %seems to be a problem with the CIP_counts - first bin
    %is not looking right - 2nd bin is what should be the 1st etc.
    sVol_CIP(i+1,:)=SampleVolumeCIP(TAS_sVol,i,1); %TAS_all is the airspeed in m/s
    %on the same time base as the CIP counts
end

%actually this is done in cas_sample_volume_and_stats2 too... called CIP_total_number
total_N_CIP=sum(CIP_counts_all(:,2:end)./sVol_CIP(:,1:end)',2);

switch tas_type
    case 'constant'
        number_factor = 60./(air_speed_1D/100)'; %multiply by this to get the corrected numbers
                                 %sample volume is proportional to the airspeed so number is inverslely proportional
%number_factor is based on CAS time (it is interpolated from CIP_time_all to CAS_time_all in 
%cas_sample_volume_and_stats

        total_N_CIP = total_N_CIP.*number_factor;
end
                                 
