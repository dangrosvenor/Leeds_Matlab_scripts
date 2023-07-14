%take the updraught field of a flight segment as a function of x
%Then assume that the updraught is constant with height and integrate backwards in time
%to get to the first position upstream to give an upstream profile.
%Will start with the last 35 km or so of leg 1 as before this we have some NaN data 
f=1e6*28.97/18; %conversion between MR and ppmv - use 18 for water vapour and 48 for ozone


%need to use CAS data to better determine when have LWC (so when saturated)
%Plus add the temperature change during the backwards changes.


angle_offset=0; %angle between the flight direction and the wind direction
%set to zero to switch this off
ivar_angle_offset=0; %flag to calculate a different offset angle for each location
%relative to a fixed bearing
offset_bearing=20; %this is the angle from north FROM which the vector points
  %like a wind direction
ivert_slope=1; %flag to use the vertical slope concept
vert_slope_angle = 68.2; %angle of vertical slope 

ionly_sat=0; %flag to ignore any points where are not using the saturated values to show where they are
iuse_sat=1; %flag to use saturated values if there is LWC or RH>1
iuse_qv_liq=1; %flag to use the flight track vapour values as calculated using
      %the 'liq' flag instead of 'ice' for the calc from the dew point
icas_lwc_presence=1; %flag to use the CAS for the determination of LWC rather than the hotwire
                    %as is less noisy (actual but values likely wrong).
iuse_cas_lwc=1; %flag to use the CAS LWC as the aircraft LWC
%for starting trajectories and for when cross the track
RH_thresh=1;
RH_thresh=0.85; %threshold for RH above which to assume saturation and use the 
%saturated qv value (from temperature) rather than the measured (prob unreliable) vapour value

%W_offset=1.5;
W_offset=0;

ioverride_U_above=0;
U_override=28;  %was 20.11 m/s orginally (=Usound(end))

Uextend_shear=0; %uses a shear instead of constant value - doesn't use the values above

iUconstant = 1; %flag to use constant horizontal velocity
Uconstant = 20;

eval(['dat_flt = dat_flt' flight_no ';']);  %put the data for the required flight here
eval(['time_flt2 = time_flt' flight_no ';']); %put time for the flight in time_flt2

%get col_alt, etc.
set_column_numbers_for_flight_data

Aircraft_position_and_trajecctory %run load_WRF_vars first

%set the time period for the bit of the the filght to work with
% xlimits=[20+20/60 20+40/60]; %the first leg of flight 102
% [inds_1 inds_2] = findheight_nearest(time_flt2,xlimits(1),xlimits(2));
% inds=inds_1:inds_2;

%calculate the distances flown by the plane
eval(['X_flt = X_flt' flight_no ';']);
eval(['Y_flt = Y_flt' flight_no ';']);
dist = [0; cumsum(sqrt( (diff(X_flt)).^2 + (diff(Y_flt)).^2 ))];
%this is the cumulative distance along the flight track

%get the data for the flight
get_flight102_envir_sounding; %get the sounding profile measured at the end of the flight
%will use this for the pressure profile and will interpolate from this based on height
%this puts the profile in Psound and the height in Zsound

%also gives a U profile - will use this instead of the along flight track values since
%they are quite variable and possibly wrong
%will smooth the profile first and extend it to higher altitudes
nfilter=5;
zmax=6000;
dz_grid=10;
zextend = [Zsound(end)+dz_grid:dz_grid:zmax];
Zsound_extend = [Zsound' zextend];

if ioverride_U_above==1
    Uend=U_override;
else
    Uend=Usound(end);
end

Uz_lower = 2800;
Uz_lower_U = 11;
Uz_upper = 3100;
Uz_upper_U = 20;
Ugrad = (Uz_upper_U-Uz_lower_U) / (Uz_upper - Uz_lower); %gradient/wind shear

if Uextend_shear==1
    Usound_extend = [Usound' Usound(end)+(zextend-zextend(1))*Ugrad];
else
    Usound_extend = [Usound' Uend*ones(size(zextend))];
end

if iUconstant==1 %constant U
    Usound_extend = [Uconstant*ones(size(Zsound_extend))];
end

[Zsound_smooth,Usound_smooth] = smooth_data(Zsound_extend,Usound_extend,nfilter);

%now construct a temperature profile for calculating the pressure profile since we
%do not have pressure data above 3200m - assume constant lapse rate
Zsound_lower_grad = 2450; %3166  %lowest height for calculating the gradients
Zsound_lower=2450;
Tsound_lower=interp1(Zsound,Tsound,Zsound_lower);    %degC;
Psound_lower=interp1(Zsound,Psound,Zsound_lower);    %hPa

Zsound_upper = 3050; %3080; %3166
Tsound_upper=interp1(Zsound,Tsound,Zsound_upper);    %-11.24;

lapse_rate = (Tsound_upper-Tsound_lower)/(Zsound_upper-Zsound_lower_grad); %degrees/metre
Zsound_pressure = [Zsound_lower:dz_grid:zmax];

T0=Tsound_lower;
Tend=(Zsound_pressure(end)-Zsound_lower)*lapse_rate + T0; %assume a constant lapse rate to extrapolate to higher alts
Tsound_pressure=273.15+[T0 : (Tend-T0)/(length(Zsound_pressure)-1) : Tend]; %constant lapse rate from Zsound_lower to zmax

%will calculate P because don't have the full range of P up to heights will need - assume hydrostatic
%But, need to assume a tempeature profile - use constant lapse rate approximation
ZSPAN=Zsound_pressure; %can specify all the points we want to get for the solution.
%perhaps consider increasing the resolution of h for consistency with the sounding
[h,P_hyd] = ODE45(@hydrostatic,ZSPAN,Psound_lower*100,[],Zsound_pressure,Tsound_pressure); %solve hydrostatic equation - uses TSPAN to interpolate temperautre for a given H
Psound_extended = interp1(h,P_hyd,Zsound_pressure); %interpolate T onto the h-grid determined by the ODE solver

%might want to consider calculating the pressure change for each parcel. We know P,T at at given
%point and given dz we can calculate dp from the hydrostatic equation assuming little change
%in P,T over the interval. The new temperature can then be calcualted from the moist or dry
%adiabtic ascent routine for a given dp. This would probably mean that pressure varied horizontally

        
        


W = w2_turb(:); % (m/s)
U = dat_flt(:,col_wind); %(m/s) will need this to calculate the movement of an air parcel through
%the streamlines
winddir = dat_flt(:,col_winddir)+180; %winddir (degrees)
winddir(winddir>360)=winddir(winddir>360)-360;
Z = dat_flt(:,col_alt); %altitude (m)
T = 273.15+dat_flt(:,col_temp);  %temperature (K)
qv = qv_flt102_humi;   %vapour (kg/kg) humicap
%qv = qv_flt102_fp;   %vapour (kg/kg) hygrometer
LWC = interp1(CIP_time_all/3600,LWC_CAS_all,time_flt2); %LWC (g/m3) - interpolate for dat_flt time base
dewT_humi = dat_flt(:,col_frostpoint_humi);
dewT_hygro = dat_flt(:,col_frostpoint_hygro);

P = 100*dat_flt(:,col_press); %pressure (Pa)


%get the CAS LWC
itimser='CAS plots';
man_choose_itimser=1;
man_choose_flt_graph=1;
i_highlight_path=1;       
time_graph = {'LWC_dist_CAS'};
instrument={'BAS CAS'};
subplotting=0;
ix_distance=1;
air_speed_type = 'CIP probe';
airspeed_constant=0;
timeseries_dan   
        
LWC_cas = interp1(CAS_time_all/3600,LWC_dist_cas,time_flt2); %LWC (g/m3) - interpolate for dat_flt time base


x_start=134.4; %km - the start of the bit of leg one that we want to work with
x_end=156.9; %km - the end of the bit of leg one that we want to work with
%note for this flight leg the plane was moving away from the mountain
%thus increasing dist values indicate we are moving downstream
[inds_1 inds_2] = findheight_nearest(dist,x_start,x_end);
inds2=inds_1:inds_2;

W2=W(inds2)+W_offset;  %+0.5 Apply offset to W here if required
U2= U(inds2); %U should be set positive if flying downstream and negative if upstream
winddir2 = winddir(inds2); %winddir
dir_flt2 = dir_flt(inds2); %direction (trjectory) of the aircraft
if ivar_angle_offset==1
    X_origin = X_flt(inds2) - X_flt(inds2(1)); %x-distance from start pos
    Y_origin = Y_flt(inds2) - Y_flt(inds2(1)); 
    %will use the start pos as a reference for the north angle
    [ilat,ilon] = getind_latlon_quick(lat2d.var,lon2d.var,dat_flt((inds2(1)),col_lat),dat_flt((inds2(1)),col_lon),0.1);
    %calculate angle between points and the loc of the first point
    for iloc=1:length(X_origin);
        dir_flt_origin(iloc)=wind_dir_compass_from_uv_wrf(X_origin(iloc),Y_origin(iloc),lat2d,lon2d,ilat,ilon,DX,DY);
    end
    angle=dir_flt_origin - offset_bearing; %angle between aircraft dir and assumed constant wind dir
    distX = sqrt(X_origin.^2+Y_origin.^2)';
else
    angle=angle_offset; %constant offset
    distX = dist(inds2);
end

X2=distX.*1000.*cos(angle*pi/180); %convert to m
Z2=Z(inds2);
T2=T(inds2);
P2=P(inds2);
qv2=qv(inds2);
qv_used=qv2;
rho=density(P2,T2);
LWC2=LWC(inds2)/1000 ./rho; %convert to kg/kg
LWC2(LWC2<0)=0;

LWC2_cas=LWC_cas(inds2)/1000 ./rho; %convert to kg/kg
LWC2_cas(LWC2_cas<0)=0;

dewT_humi2 = dewT_humi(inds2)+273.15;
dewT_hygro2 = dewT_hygro(inds2)+273.15;
qv_liq_humi2 = SatVapPress(dewT_humi2,'goff','liq',P2,1)/f;
qv_liq_hygro2 = SatVapPress(dewT_hygro2,'goff','liq',P2,1)/f;

qsat_flight = SatVapPress(T2,'goff','liq',P2,1)/f;
equiv_flight = equivalent_potemp(T2,P2,qsat_flight);

dz_max=5;  %max vertical distance (m) to move in dt
LWC_thresh=0.057/1000; %kg/kg

if icas_lwc_presence==0
    LWC_thresh=0.1/1000; %kg/kg
else
    LWC_thresh=0.01/1000; %kg/kg  - can use a much lower threshold as the CAS
    %is much less noisy and therefore better at determining the presence of LWC
end

%do the tilted gwave field calculation
calculate_uv_field_vertical_angle_gwaves

dt = - dz_max / max(W2); %make negative for running backwards

clear zz tot2 equiv2 temp2 L2
ip2=0;
for ip=1:1:length(W2)  %loop over all the points observed by the aircraft - will do the
                     %backwards intergration over all these points
                     
     ip2=ip2+1; %counter for points that have integrated for
               
     fprintf(1,'ip=%d/%d\n',ip,length(W2));
     %initialise starting positions for the datapoint   
     x = X2(ip); % (metres)
     z = Z2(ip);
     %     p = P2(ip); %pressure Pa
     p = interp1(Zsound_pressure,Psound_extended,z); %use a profile for pressure instead of the 
     %observed value at the moment

     
     temp = T2(ip);
     
     if iuse_qv_liq==1
         q = qv_liq_humi2(ip);
     else
         q = qv2(ip); %kg/kg
     end
     
     if iuse_cas_lwc==1
         L=LWC2_cas(ip);
     else
         L = LWC2(ip); %kg/kg
     end
     
     if icas_lwc_presence==1
         L_present = LWC2_cas(ip); %LWC from the CAS - probably a better indicator
     %of the presence of liquid (less noisy) than the hotwire even if the value
     %may not be correct
     else
         L_present = L;
     end
     
     qsat=SatVapPress(temp,'goff','liq',p,1)/f;
     RH = q/qsat;
     
     %%% *** change to take into account CAS
     if L_present>LWC_thresh | RH>RH_thresh %if have liquid water then use the saturated MR value as is likely to be more accurate         
         if iuse_sat==1
             q=qsat;
         end
         qv_used(ip)=qsat; %save for record of the flight track values
         qv_used_flag(ip)=1;
     else
         qv_used_flag(ip)=0;
          if ionly_sat==1
             q=NaN; %test to see what things look like for just the saturated points
                     %i.e. if don't rely on the vapour instrument
          end
     end
     tot = q+L; %total water
     
     equiv = equivalent_potemp(temp,p,q);
     
     tot_X2(ip) = tot;
     

     t=0; %time
     
     x_target = X2(1); %keep going until reach the most upstream position of the segment
     %since the plane was moving downstream this is the first point
     %now start the backwards integration in time
     i=1;
     updated=0;
     sign_old=1;
     stop_while=0;
     while stop_while==0 & x(i)>x_target
         t=t+dt;
         i=i+1;
         
         if isnan(x(i-1))==1 | isnan(z(i-1))==1
             stop_while=1;
             continue
         end
             
         
         if ivert_slope==1 %flag to use the vertical slope concept
             u(i-1) = interp2(X2,heights_map,U_af',x(i-1),z(i-1)); %use a value as a function of height only
%             u(i-1) = 15;
             w(i-1) = interp2(X2,heights_map,W_af',x(i-1),z(i-1)); 
         else             
             u(i-1) = interp1(Zsound_extend,Usound_extend,z(i-1)); %use a value as a function of height only
             w(i-1) = interp1(X2,W2,x(i-1));
         end
         
         x(i) = x(i-1) + u(i-1)*dt; %perhaps consider using contant U as get some strange variations?
         %also consider using the wind direction to calculate a component in a given direction?
         z(i) = z(i-1) + w(i-1)*dt;
         

         
         %interpolate pressure based on altitude and the current height
         p(i) = interp1(Zsound_pressure,Psound_extended,z(i));
         
         tot(i)=tot(i-1); %conserve total water
         
         qsat = SatVapPress(temp(i-1),'goff','liq',p(i-1),1)/f;
         
         %calculate temperature change depending on whether are in cloud or not
         rho=density(p(i-1),temp(i-1));
         if L(i-1)>1e-12 | q(i-1) >= qsat
             [LWC_ad,T_ad]=adLWC_PaulLawson_simple(p(i-1),temp(i-1),[p(i-1)-1 p(i)]); %move from p(i-1) to p(i)
             LWC_ad_kg = LWC_ad(2)/1000/rho; %convert from g/m3 to kg/kg
             L(i) = L(i-1)+LWC_ad_kg; %add to old LWC
             if L(i)<0 %when have evaporated all the LWC avaiable
                 L(i)=0;
                 q(i)=q(i-1); %same q as before
             else
                 q(i) = q(i-1)-LWC_ad_kg;  %remove LWC formed - LWC will be negative
                            %for evaporation - but only if had LWC left!
             end
                          
             temp(i) = T_ad(2) + 273.15;
         else
             q(i) = q(i-1); %no change in qv
             L(i) = L(i-1); %should already by >=0
             th = temp(i-1) * (1e5/p(i-1)).^0.286; %Dry adiabat
             temp(i) = th / (1e5/p(i)).^0.286;
         end
                           
         zplane = interp1(X2,Z2,x(i)); %where the aircraft was for this x
         sign_new = sign(z(i)-zplane);
         if (sign_new~=sign_old & i~=2);%if have crossed the aircraft path then change to the value observed at that point
             %observed value at the moment


             temp(i)=interp1(X2,T2,x(i));     
             
             if iuse_qv_liq==0
                 q(i) = interp1(X2,qv2,x(i));
             else
                 q(i) = interp1(X2,qv_liq_humi2,x(i));
             end
             
             if iuse_cas_lwc==1
                 L(i)= interp1(X2,LWC2_cas,x(i));
             else
                 L(i)=interp1(X2,LWC2,x(i));
             end
             
             if icas_lwc_presence==1
                 L_present = interp1(X2,LWC2_cas,x(i));
             else
                 L_present = L(i);
             end
             
             qsat=SatVapPress(temp(i),'goff','liq',p(i),1)/f;
             RH = q(i)/qsat;
%             if L(i)>LWC_thresh | RH>1 %if have liquid water then use the saturated MR value as is likely to be more accurate
             if L_present>LWC_thresh | RH>RH_thresh %if have liquid water then use the saturated MR value as is likely to be more accurate
                  if iuse_sat==1
                     q(i)=qsat;
                  end
             else
                 if ionly_sat==1
%                     q(i)=NaN; %test to see what things look like for just the saturated points
                 end
                 %i.e. if don't rely on the vapour instrument
             end
             tot(i) = q(i)+L(i); %total water
             
             updated=1;

         end

         sign_old=sign_new;   
         equiv(i) = equivalent_potemp(temp(i),p(i),q(i));
         
     end %end of while loop for a particular starting point
     
     inot_nan = find(isnan(x)==0);
     %NEED to clear all of these
     if length(x(inot_nan))>1 
         zz(ip2,:) = interp1(x(inot_nan),z(inot_nan),X2); %interpolate onto the regular X2 grid     
         tot2(ip2,:) = interp1(x(inot_nan),tot(inot_nan),X2);
         equiv2(ip2,:) = interp1(x(inot_nan),equiv(inot_nan),X2);      
         temp2(ip2,:) = interp1(x(inot_nan),temp(inot_nan),X2);       
         L2(ip2,:) = interp1(x(inot_nan),L(inot_nan),X2);                
     else
         zz(ip2,:) = [z(1) ones([1 length(X2)-1])*NaN];
%         zz(ip2,2:length(X2))=NaN;
         tot2(ip2,:) = [tot(1) ones([1 length(X2)-1])*NaN];
         equiv2(ip2,:) = [equiv(1) ones([1 length(X2)-1])*NaN];         
         temp2(ip2,:) = [temp(1) ones([1 length(X2)-1])*NaN];                  
         L2(ip2,:) = [L(1) ones([1 length(X2)-1])*NaN];                           
     end
     
end  %end of loop over all starting points

clear tot3 equiv3 temp3 L3
%interpolate onto a regular z grid
for ix=1:length(X2)
    inan=isnan(zz(:,ix));
    inani=find(inan==0);
    if length(inani)>1
        tot3(ix,:)=interp1(zz(inani,ix),tot2(inani,ix),Zsound_pressure);
        equiv3(ix,:)=interp1(zz(inani,ix),equiv2(inani,ix),Zsound_pressure);    
        temp3(ix,:)=interp1(zz(inani,ix),temp2(inani,ix),Zsound_pressure);            
        L3(ix,:)=interp1(zz(inani,ix),L2(inani,ix),Zsound_pressure);                    
    else
        tot3(ix,:)=ones(size(Zsound_pressure))*NaN;
        equiv3(ix,:)=ones(size(Zsound_pressure))*NaN;        
        temp3(ix,:)=ones(size(Zsound_pressure))*NaN;    
        L3(ix,:)=ones(size(Zsound_pressure))*NaN;            
    end
end





