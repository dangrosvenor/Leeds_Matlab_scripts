filename_acpim_sounding = 'C:\cygwin\home\dan\acpim1\acpim\TEXT\wave_cloud_case_sounding_BACK_TRAJ_1.1.txt';
no_write=0; %don't write out the file - just calcualte the fields

write_velocity_file=1;
filename_acpim_velocity = 'C:\cygwin\home\dan\acpim1\acpim\TEXT\vertical_vel_plus_1.5.txt';


%%% flag to extrapolate the sounding to reach the surface as a test in ACPIM
%continues an idealised lapse rate to the ground and assumes constant RH
%means that need to decrease the temperature towards the ground in order
%to get 1e5 hPa there.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iextrapolate_to_surface=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iz0=0;  %set the z-axis to start at z=0, i.e. relative to the bottom of the sounding
        %Will be keeping the orginial sounding, but just changing z - need to change
        %in Dynamics.c too.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


z2_sound=z;
T_h_sound=T_h; %temperature (K)
qv2_sound=qv_prof_humi;
P2_sound=P_wave;

%apply the prescribed constant qv to the heights beyond the top of the sounding
%if this was done in Twave
if iuse_mean_extrap==1;
    ih=findheight(z2_sound,max(Zsound));
    qv2_sound(ih+1:end) = qvhumi_mean;
end

iuse_constant_RH_abv_lift=1;

if iuse_constant_RH_abv_lift==1
%    ih=findheight(z2_sound,z0damp+dz_damp);
    ih=findheight(z2_sound,max(Zsound));     %here am setting this to start at the top of the sounding
                                       %as otherwise the vapour profile looks a bit silly
    qsat_humi = satvappress(T_h_sound(ih+1:end),'goff','liq',P_wave(ih+1:end),1)/f;
    RH = qv2_sound(ih+1)'./qsat_humi(1); %constant RH to apply
    
    RH=0.2; %override
    
    qv2_sound(ih+1:end) = RH * qsat_humi';
end





if iextrapolate_to_surface==1
    %find a lapse rate - was fairly constant over these heights
    Z1=895;
    z2_sound=1100;
    iZ=findheight(Zsound,Z1);
    iz2_sound=findheight(Zsound,z2_sound);
    lapse_rate=(Tsound(iz2_sound)-Tsound(iZ) ) / (z2_sound-Z1);
    lapse_rate=lapse_rate*-4;

    %apply this lapse rate to the part below where we already have
    dz=min(diff(z));
    z2_sound=[0:dz:max(z)];
    iz=findheight(z2_sound,z(1));
    T_h_sound=[((z2_sound(1:iz-1) - z2_sound) * lapse_rate + Tsound(iz2_sound)+273.15 )' ; T_h];
    
    %now solve the hydrostatic equation from this point down to the surface
    [Hhyd,Phyd] = ODE45(@hydrostatic,[z2_sound(iz) z2_sound(1)],P_wave(1),[],z2_sound,T_h_sound);
    
    %need to interpolate onto the z2_sound grid as ODE45 creates its own grid.
    P2_sound = [interp1(Hhyd,Phyd,z2_sound(1:iz-1))'; P_wave];
    
    %keep RH constant
    RH = qv_prof_humi(1)/( SatVapPress(T_h(1),'goff','liq',P_wave(1),1)/f);
    qv2_sound = [RH*( SatVapPress(T_h_sound(1:iz-1),'goff','liq',P2_sound(1:iz-1),1)/f); qv_prof_humi'];
    
end
    
no_points=length(T_h);

if iz0==1
    z2_sound=z2_sound-min(z2_sound);
end

if no_write==1
            theta = T_h_sound .* (1000e2./P2_sound).^0.286;
else
    
fid=fopen(filename_acpim_sounding,'wt');

fprintf(fid,'%d\n',no_points);

sound_case =1;
switch sound_case
    case 0 %for this type ACPIM reads in pressure, temperature, dew point and height
        dew=Tdew(qv2_sound,P2_sound)-273,15;
        for i=1:no_points
            fprintf(fid,'%f %f %f %f\n',P2_sound(i)/100,T_h_sound(i)-273.15,dew(i),z2_sound(i));
        end
        
    case 1  %for this type ACPIM reads in height, potemp and vap MR
        theta = T_h_sound .* (1000e2./P2_sound).^0.286;
        fprintf(fid,'%s\n','scratch'); %needs this as reads it in       
        for i=1:no_points
            fprintf(fid,'%f %f %f %f\n',z2_sound(i),P2_sound(i),theta(i),qv2_sound(i));
        end
end

fclose(fid);

end


if write_velocity_file==1
    fid=fopen(filename_acpim_velocity,'wt');
    no_points=length(vert_vel);
    fprintf(fid,'%d\n',no_points);
    
    for i=1:no_points
            fprintf(fid,'%f %f\n',time_wave(i),vert_vel(i));
    end
        
end

disp('Done write sounding');

