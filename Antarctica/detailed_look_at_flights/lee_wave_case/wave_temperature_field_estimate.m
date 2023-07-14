%from aircraft profile 21:40:38 to 22:01:41 (at the end of the flight)
%this defines the lapse rate above the inversion at around 2300-2400 m
%lapse rate above here is approx constant
ijust_calculate_for_sounding=1;  %flag to say that we just want to compute values for the write of the
%sounding and to not do the time consuming interpolation onto regular grids, plotting, etc.
iplot=0;    
filedir='C:\Documents and Settings\dan\My Documents\logbook\Antarctica\Flights and instruments_Feb2010\Lee wave cloud case\Wave_cloud_Matlab_plots\';
filedir='Y:\BAS_flights\flight102\Wave_model_ACPIM_plots\Matlab_model\';

isave=0;

warning off

iprofs=1;
vap_factor=1;
%vap_factor=1.1; %factor to increase the vapour by for testing (to get a cloud
%base more consistent with that required to give the hotwire LWCs).
add_constant_temp=-0.1;
add_constant_eq=0.5;

add_constant_temp=0;
add_constant_eq=-0.25;

%add_constant_temp=0;
%add_constant_eq=0;

%value to add to the top of the sounding for waves starting higher up (unknown from
%the measurements)
iuse_mean_extrap=0; %whether to extrapolate qv to alitudes we don't have data for by assuming a
            %constant value
qvfp_mean = 2.2e-3; %1.1e-3; %mean(qfp_sound(isound_lower:isound_upper))/2;
qvhumi_mean = 2.2e-3; %1.1e-3; %mean(qhumi_sound(isound_lower:isound_upper))/2;

iuse_constant_RH_above_sounding=1; %NOTE - this overrides the above qvfp_mean
                                   %and qvhumi_mean values with a constant
                                   %specified RH 



iuse_equiv_constant_lapse_rate=1; %flag to say that want to assume a constant lapse rate for equiv potemp
%(and for temperature) - then also calculates the water vapour from the temperature (constant lapse rate)
% and pressure (from hydrostatic equation). Water vapour at the top of the profile will be replaced by
% the applied upper value (qvfp_mean and qvhumi_mean)

iuse_constant_temp_lapse_rate=1; %flag to say whether to assume a constant temeprature lapse rate for the part
%of the profile covered by the sounding. In all cases the profile is extended above there using a constant lapse
%rate (for the pressure calculation?)

svp_method_liq='goff'; %previously used
svp_method_liq='buck2'; %changed for consistency with ACPIM
svp_method_ice='goff'; %previously used
svp_method_ice='murphy'; %changed for consistency with ACPIM


    
if vap_factor>1
    fprintf(1,'\n**** WARNING! Applying a factor of %f to the vapour sounding!! ****\n',vap_factor);
end
if add_constant_temp>0
    fprintf(1,'\n**** WARNING! Adding a temperature of %f to the lapse rates!! ****\n',add_constant_temp);
end
if add_constant_eq>0
    fprintf(1,'\n**** WARNING! Adding a temperature of %f to the lapse rates!! ****\n',add_constant_eq);
end



z0=3000; %height of the zero displacement line of the largest wave
zmax=6000;  %max height of the z-grid


f=1e6*28.97/18; %conversion between MR and ppmv - use 18 for water vapour and 48 for ozone

%%%% ****  get the environmental (ascent & descent) sounding ****
get_flight102_envir_sounding
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%will only consider the approx constant lapse rate of the upper part of this sounding (above 2420m)
% isound_lower=findheight_nearest(Zsound,2420);
% Zsound_lower = Zsound(isound_lower); %2422
% Tsound_lower=Tsound(isound_lower);    %-7.8;
% Psound_lower=Psound(isound_lower);    %-7.8;
% eq_sound_lower_humi=equiv_sound_humi(isound_lower);    %
% eq_sound_lower_fp=equiv_sound_fp(isound_lower);    %


Zsound_lower_grad = 2450; %3166  %lowest height for calculating the gradients
Zsound_lower=2450;
Tsound_lower=interp1(Zsound,Tsound,Zsound_lower);    %-11.24;
eq_sound_lower_humi=interp1(Zsound,equiv_sound_humi,Zsound_lower);    %297.16
eq_sound_lower_fp=interp1(Zsound,equiv_sound_fp,Zsound_lower);    %297.11

override_grad=0;
if override_grad==1
    eq_sound_lower_humi = 297.16+0.5;    %289.6
    eq_sound_lower_fp = 297.11+0.5; 
    Zsound_lower=2500;
end
        
Psound_lower=interp1(Zsound,Psound,Zsound_lower);    %


%[maxZ isound_upper]=max(Zsound); 

%Zsound_upper = Zsound(isound_upper); %3166
%Tsound_upper=Tsound(isound_upper);    %-11.24;
%eq_sound_upper_humi=equiv_sound_humi(isound_upper);    %
%eq_sound_upper_fp=equiv_sound_fp(isound_upper);    %

Zsound_upper = 3050; %3080; %3166
Tsound_upper=interp1(Zsound,Tsound,Zsound_upper);    %-11.24;
eq_sound_upper_humi=interp1(Zsound,equiv_sound_humi,Zsound_upper);    %
eq_sound_upper_fp=interp1(Zsound,equiv_sound_fp,Zsound_upper);    %


ilapse_rate = 'ascent';
ilapse_rate = 'descent';

switch ilapse_rate
    case 'descent'
        lapse_rate = (Tsound_upper-Tsound_lower)/(Zsound_upper-Zsound_lower_grad); %degrees/metre
        lapse_rate_eq_humi = (eq_sound_upper_humi-eq_sound_lower_humi)/(Zsound_upper-Zsound_lower_grad); %degrees/metre
        lapse_rate_eq_fp = (eq_sound_upper_fp-eq_sound_lower_fp)/(Zsound_upper-Zsound_lower_grad); %degrees/metre        
    case 'ascent'
        Tsound_lower = -5.5;
        Zsound_lower = 2450;
        lapse_rate = ( -11.4 - Tsound_lower ) / ( 3147 - Zsound_lower) ;
        Psound_lower = 722.3;
end
%once the lapse rate is calulated we don't use Zsound, etc anymore.

i_use_ascent_for_qv=0;
if i_use_ascent_for_qv==1
    %qfp_sound = interp1(Zsound2,qfp_sound2,Zsound);
    %qhumi_sound = interp1(Zsound2,qhumi_sound2,Zsound);
    qfp_sound = qfp_sound2; 
    qhumi_sound = qhumi_sound2;
    Zsound = Zsound2;  %Zsound is used in Twave (only for interpolation of the vapour profile)
end






% zsound_upper = 3166;
% Tsound_upper=-11.24;
% Psound_upper=659.2*100; %Pa
% qfp_sound_upper=2.3589*1e-3; %qv from frost point hygrometer (kg/kg)
% qhum_sound_upper=2.4095*1e-3; %qv from frost point hygrometer (kg/kg)

%wave crests are at 3100-3300 m.
%from the end of flight profile the air above the inversion at ~2400m is much moister than that below
%thus this could explain the presence of the waves here - might indicate that the wave field could be uniform
%and that it is the presence of this moist air that forms clouds here. Although both the presence of the
%waves and the moist layer may be related to the inversion in which case may still have a damped wave field

%BUT - Tom mentioned that the hygrometers have a problem when close to 100% RH so the jump could be due to that?




formulation='new'
switch formulation
    case 'new'
        dz_grid=10;
        z = [Zsound_lower:dz_grid:zmax];
        
        T0=Tsound_lower;
        Tend=(z(end)-Zsound_lower)*lapse_rate + T0; %assume a constant lapse rate to extrapolate to higher alts
        T_wave=273.15+[T0 : (Tend-T0)/(length(z)-1) : Tend] + add_constant_temp;
        
        if iuse_constant_temp_lapse_rate==0
%            iz_0=findheight(z,Zsound_lower); (=1)
            iz_1=findheight(z,Zsound_upper);            
            T_wave(1:iz_1)=interp1(Zsound,Tsound+273.15,z(1:iz_1));
        end
        
        eq_0=eq_sound_lower_humi;
        eq_end=(z(end)-Zsound_lower)*lapse_rate_eq_humi + eq_0; %assume a constant lapse rate to extrapolate to higher alts                
        eq_humi=[eq_0 : (eq_end-eq_0)/(length(z)-1) : eq_end] + add_constant_eq;        
        
        eq_0=eq_sound_lower_fp;
        eq_end=(z(end)-Zsound_lower)*lapse_rate_eq_fp + eq_0; %assume a constant lapse rate to extrapolate to higher alts                
        eq_fp=[eq_0 : (eq_end-eq_0)/(length(z)-1) : eq_end] + add_constant_eq;        

        iextend_lapse_below=1;
        lowest_height=2400; %max height for this is Zsound_lower-dz_grid
        if iextend_lapse_below==1
            z_below = [lowest_height:dz_grid:z(1)-dz_grid];
            %ilowest = findheight_nearest(z,lowest_height);
            T_wave = [interp1(z,T_wave,z_below,'linear','extrap') T_wave];
            eq_humi = [interp1(z,eq_humi,z_below,'linear','extrap') eq_humi];            
            eq_fp = [interp1(z,eq_fp,z_below,'linear','extrap') eq_fp];                        
            press0 = 100*interp1(Zsound,Psound,z_below(1));
            z=[z_below z];
            Zsound_lower=z(1);
        else        
            press0=Psound_lower*100;
        end
        
        iextend_below=1;
        lowest_height=1000;
        if iextend_below==1
            z_below = [lowest_height:dz_grid:z(1)-dz_grid];
            %ilowest = findheight_nearest(z,lowest_height);
            T_wave = [interp1(Zsound,273.15+Tsound,z_below) T_wave];
            eq_humi = [interp1(Zsound,equiv_sound_humi,z_below) eq_humi];            
            eq_fp = [interp1(Zsound,equiv_sound_fp,z_below) eq_fp];                        
            press0 = 100*interp1(Zsound,Psound,z_below(1));
            z=[z_below z];
            Zsound_lower=z(1);
        else        
            press0=Psound_lower*100;
        end
        
        %will calculate P because don't have the full range of P up to z(end) - assume hydrostatic
%        ZSPAN=[z(1) z(end)];
        ZSPAN=z; %can specify all the points we want to get for the solution.
        %perhaps consider increasing the resolution of h for consistency with the sounding
        [h,P_wave] = ODE45(@hydrostatic,ZSPAN,press0,[],z,T_wave); %solve hydrostatic equation - uses TSPAN to interpolate temperautre for a given H
        T_h = interp1(z,T_wave,h); %interpolate T onto the h-grid determined by the ODE solver
                
        pot=T_h.*(1000e2./P_wave).^0.286;
        N=sqrt( 9.81/pot(1) * (pot(end)-pot(1))/(h(end)-h(1)) )
        


            
         
    case 'old'
        z=[-1:0.01:1]*1000;
        z_actual = z_base + z - z(1);
        dZ=z-z(1);
        pot0=295;
        %assume constant lapse rate temperature starting at pot0
        T0=-11;
        Tend=-18;
        %Tend=-25;
        T_wave=273.15+[T0 : (Tend-T0)/(length(z)-1) : Tend];
        press0=((273+T0)/pot0).^(1/0.286) * 1000e2;       

        ZSPAN=[z(1) z(end)];
        [h,P_wave] = ODE45(@hydrostatic,ZSPAN,press0,[],z,T_wave); %solve hydrostatic equation - uses TSPAN to interpolate temperautre for a given H
        T_h = interp1(z,T_wave,h);

        pot=T_h.*(1000e2./P_wave).^0.286;
        %N=0.01; %assume constant stratifcation
        N=sqrt( 9.81/pot(1) * (pot(end)-pot(1))/(h(end)-h(1)) )

   

end
h_wave_model = h;

%interpolate equivalent potemp onto the h vertical grid


if iuse_equiv_constant_lapse_rate==1
    eq_h_humi = interp1(z,eq_humi,h);
    eq_h_fp = interp1(z,eq_fp,h);
    qv_prof_humi = find_qv_from_equiv_potemp_P_and_T(eq_h_humi,P_wave,T_h);
    qv_prof_fp = find_qv_from_equiv_potemp_P_and_T(eq_h_fp,P_wave,T_h);    
else %default option
    qv_prof_humi = interp1(Zsound,qhumi_sound,h);
    inan=find(isnan(qv_prof_humi)==1);
%    qv_prof_humi(inan(1):end)=qvhumi_mean;
    
    qv_prof_fp = interp1(Zsound,qfp_sound,h);
%    qv_prof_fp(inan(1):end)=qvfp_mean;
    eq_h_humi = equivalent_potemp(T_h,P_wave,qv_prof_humi);
    eq_h_fp = equivalent_potemp(T_h,P_wave,qv_prof_fp);

end



x=[0:0.01:9]; %km
xf=20/9 * (1/3) ;  %presumably in km^{-1}
%xf=20/9 * (1/3) * 0.8 ;  %presumably in km^{-1}
xf=2*pi/7.5;  %y=1000*yf*sin(kx), with x in km. For x=wavelength sin(kx)=sin(0) so k*wavelength=2pi
%then xf=k=2pi/wavelength. And wavelength was around 7.5 km.

yf=1/4 *2.2;
yf=1/4 *1.1;  %2*yf is the peak to peak amplitude
yf=1/4 *1.5;
%yf=1/4 *2.3;
yf=0.623;
%yf=0.7;
%yf=0.772; %km  : wmax=11.43 m/s 
yf=0.405; %km  : wmax=6 m/s
%yf=0.6075; %km : wmax=9 m/s
%yf=0.55;
%yf=0.8;
%yf=0.744; %based on max displacment from the plot (estimated from the equiv
%potemp plot) - corresponds to wmax = 11.0 m/s. Parcel observed at 3.229 and
%so started at 2.455 km

WS=20; %wind speed in m/s since dy/dx is in km/km and so dimensionless

wave_type = 'sin';
%wave_type = 'asymmetrical';
switch wave_type
    case 'sin'
        y=1000*sin(x*xf)*yf; %convert to m
        %dy/dx = xf*yf*cos(x*xf) %km/km
        %dy/dt = dy/dx * dx/dt = u*dy/dx
    case 'asymmetrical'
        Yoff=0.3;
        Xoff=asin(-Yoff);
        y=1000*yf* ( Yoff + sin(x*xf + Xoff) );
        
end



V = WS*yf*xf*cos(x*xf);

time=1000*x/WS; %calculate the time variation of the displacement
%if we were following the wave at the windspeed (seconds)
%doing this to calulate vertical velocities

wave_production_method='manually choose waves';
wave_production_method='automatically generate a fine mesh of waves';
switch wave_production_method
    case 'automatically generate a fine mesh of waves'

        clear TT YY VV RH_fp_ALL RH_humi_ALL RH_fp_i_ALL RH_humi_i_ALL...
            LWC_fp_ALL Tad_fp_ALL LWC_humi_ALL Tad_humi_ALL qv_wave_fp qv_wave_humi...
            PY

        nwaves=100;
        z_start=Zsound_lower*1.05; %start slightly above the lowest point we have data for
%        z_end=Zsound_upper*0.95;
         z_end = 4000;

        dz_waves = (z_end-z_start)/(nwaves-1);
        z0_waves = [z_start:dz_waves:z_end];
        
        %create a damping factor that is 1 at z0damp (i.e. no damping) and that reduces
        %to zero over a distance of dz_damp above and below
        idamp=0;
%        idamp=1;
        if idamp==1
            z0damp = 2850;
            z0damp = 2500;
%            z0damp = 2700;     
%            z0damp = 2500;  
                        
%            z0damp = 2428;  %usual choice
            
%            z0damp = 2643;

%            z0damp=3155;
%            z0damp=2400;
            
            icalc_from_m2=0;
            if icalc_from_m2==1
                h_max_neg_obs = 2615;
                disp_max_neg = 540;
                h_max_neg_z0 = h_max_neg_obs + disp_max_neg; %zero displacement height of the
                %max negative displacement observed (assuming we observed the max)
                h_neg_zero = 3165;

                m2=(disp_max_neg - 0)/(h_neg_zero-h_max_neg_z0);
                Ymax=1200;
                z0damp = h_neg_zero - Ymax/m2;            
            end

            icalc_from_m=0;
            if icalc_from_m==1
                m=(192.8+538.9)/(3.354e3-2.616e3);
                m=1.4;
                yf*1000*(m-2);
                dz_damp = -yf*1000*(m-2);

                z0damp = 3200-dz_damp;    
            end
        
 
            
            
        else
            z0damp = z0_waves; %this will set damp_factor to 1 everywhere (no damping)
        end
%        dz_damp = 900;
%        dz_damp = 1050;
        dz_damp = 1150;  %usual one
%        dz_damp = 1350;
%        dz_damp = 1250;        
%        dz_damp = 308;  
%        dz_damp = 900;
%       dz_damp = 300;
%        dz_damp = 713.5;
%        dz_damp = 100; %for severe damping
        
       
        
%        dz_damp = h_neg_zero-z0damp;

        dz_damp_exp = 500;

damp_type = 'linear';
%damp_type = 'linear only above';
%damp_type = 'linear only below';
%damp_type = 'exponential';

if idamp==0
    damp_type = 'linear';
end

switch damp_type
    case 'linear'
        damp_factor = 1 - abs(z0_waves - z0damp)/dz_damp;                
        damp_factor(damp_factor<0)=0;
    case 'exponential'
        %now dz_damp_exp becomes the height where damp_factor = 1/e
%        damp_factor = exp(-abs(z0_waves - z0damp)/dz_damp_exp);

%squared exponential
        damp_factor = exp(-((z0_waves - z0damp)/dz_damp_exp).^2);
        
    case 'linear only above'
        damp_factor = 1 - (z0_waves - z0damp)/dz_damp;                
        damp_factor(damp_factor<0)=0;
        damp_factor(damp_factor>1)=1;        
    case 'linear only below'
        damp_factor = 1 + (z0_waves - z0damp)/dz_damp;                
        damp_factor(damp_factor<0)=0;
        damp_factor(damp_factor>1)=1;            
end
        
        figure
        
        xx=damp_factor*yf*1000;
        yy=damp_factor*yf*1000+z0_waves;
        plot(damp_factor*yf*1000,damp_factor*yf*1000+z0_waves); grid
        interp1(yy,xx,3350)   %458
        close
        
        qsat_sound_top = satvappress(Tsound(end)+273.15,'goff','liq',100*Psound(end),1)/f;
        RH_constant = qhumi_sound(end) / qsat_sound_top;
        RH_constant=0.2;
        
        zmax=z0damp+dz_damp; %if idamp=0 z0damp=z0_waves and so will be a vector
        %which is why zmax(1) is used below
        for i=1:length(z0_waves)
            YY(i,:)=y*damp_factor(i)+z0_waves(i);  %wave at the correct height
            VV(i,:)=V*damp_factor(i);
            [TT(i,:),RH_fp_ALL(i,:),RH_humi_ALL(i,:),RH_fp_i_ALL(i,:),RH_humi_i_ALL(i,:)...
                ,LWC_fp_ALL(i,:),Tad_fp_ALL(i,:),LWC_humi_ALL(i,:),Tad_humi_ALL(i,:)...
                ,qv_wave_fp(i,:),qv_wave_humi(i,:),PY(i,:)]...
                =Twave(YY(i,:),h,pot,P_wave,Zsound,vap_factor*qv_prof_fp,vap_factor*qv_prof_humi,qvfp_mean,qvhumi_mean,iuse_mean_extrap...
                ,iuse_constant_RH_above_sounding,RH_constant,svp_method_liq,svp_method_ice);
            
            
            max_hdamp(i) = YY(i,1) + max(YY(i,:)) - z0_waves(i) - zmax(1);
            if z0_waves(i)>z0damp+dz_damp
                max_hdamp(i)=NaN;
            end
        end
        
        pot_wave = (Tad_humi_ALL+273.15).*(1e5./PY).^0.286;
        
        z_target = YY(65,100);
        x_target = x(100);
        fdamp = 1000*sin(x_target*xf)*yf / dz_damp;
        z0x = ( z_target - fdamp*(dz_damp + z0damp) ) / (1 - fdamp);



    case 'manually choose waves'

        [Ty,RHfp,RHhumi,RHfp_i,RHhumi_i,LWC_fp,Tad_fp,LWC_humi,Tad_humi]...
            =Twave(y+z0,h,pot,P_wave,Zsound,qv_prof_fp,qv_prof_humi,qvfp_mean,qvhumi_mean,iuse_mean_extrap);

        %can create different sin waves and then interpolate onto a regular grid

        yf2=0.5;
        yf2=1;
        offset_upper=650; %150
        offset_upper=650; %150

        y_upper=yf2*y+offset_upper;
        [Ty_upper,RHfp_upper,RHhumi_upper,RHfp_i_upper,RHhumi_i_upper,LWC_fp_upper,Tad_fp_upper,LWC_humi_upper,Tad_humi_upper]...
            =Twave(y_upper+z0,h,pot,P_wave,Zsound,qv_prof_fp,qv_prof_humi,qvfp_mean,qvhumi_mean,iuse_mean_extrap);
        V_upper = WS*yf2*(1000*yf*xf*cos(x*xf));

        %y_lower=yf2*y-offset_upper;
        y_lower = yf2*y - (z0-Zsound_lower)*0.95; %placing the start of this wave
        %as far down as possible given the data we have.

        %now work out temp at each displacement in y
        [Ty_lower,RHfp_lower,RHhumi_lower,RHfp_i_lower,RHhumi_i_lower,LWC_fp_lower,Tad_fp_lower,LWC_humi_lower,Tad_humi_lower]...
            =Twave(y_lower+z0,h,pot,P_wave,Zsound,qv_prof_fp,qv_prof_humi,qvfp_mean,qvhumi_mean,iuse_mean_extrap);
        V_lower = WS*yf2*(1000*yf*xf*cos(x*xf));

        yf2=0.35;
        yf2=1;
        offset_upper2=470; %200
        offset_upper2=870; %200
        y_upper2=yf2*y+offset_upper2;
        [Ty_upper2,RHfp_upper2,RHhumi_upper2,RHfp_i_upper2,RHhumi_i_upper2,LWC_fp_upper2,Tad_fp_upper2,LWC_humi_upper2,Tad_humi_upper2]...
            =Twave(y_upper2+z0,h,pot,P_wave,Zsound,qv_prof_fp,qv_prof_humi,qvfp_mean,qvhumi_mean,iuse_mean_extrap);
        V_upper2 = WS*yf2*(1000*yf*xf*cos(x*xf));

        yf2=0.3;
        yf2=1;
        offset_upper3=550; %220
        offset_upper3=950; %220
        y_upper3=yf2*y+offset_upper3;
        [Ty_upper3,RHfp_upper3,RHhumi_upper3,RHfp_i_upper3,RHhumi_i_upper3,LWC_fp_upper3,Tad_fp_upper3,LWC_humi_upper3,Tad_humi_upper3]...
            =Twave(y_upper3+z0,h,pot,P_wave,Zsound,qv_prof_fp,qv_prof_humi,qvfp_mean,qvhumi_mean,iuse_mean_extrap);
        V_upper3 = WS*yf2*(1000*yf*xf*cos(x*xf));


        YY=[y_upper3; y_upper2; y_upper; y; y_lower] + z0;
        TT=[Ty_upper3; Ty_upper2; Ty_upper; Ty; Ty_lower];
        VV=[V_upper3; V_upper2; V_upper; V; V_lower];
        RH_fp_ALL=[RHfp_upper3; RHfp_upper2; RHfp_upper; RHfp; RHfp_lower];
        RH_humi_ALL=[RHhumi_upper3; RHhumi_upper2; RHhumi_upper; RHhumi; RHhumi_lower];
        RH_fp_i_ALL=[RHfp_i_upper3; RHfp_i_upper2; RHfp_i_upper; RHfp_i; RHfp_i_lower];
        RH_humi_i_ALL=[RHhumi_i_upper3; RHhumi_i_upper2; RHhumi_i_upper; RHhumi_i; RHhumi_i_lower];
        LWC_fp_ALL=[LWC_fp_upper3; LWC_fp_upper2; LWC_fp_upper; LWC_fp; LWC_fp_lower];
        LWC_humi_ALL=[LWC_humi_upper3; LWC_humi_upper2; LWC_humi_upper; LWC_humi; LWC_humi_lower];

end


if ijust_calculate_for_sounding==1
    disp('Exiting wave script as ijust_calculate_for_sounding==1');
    break
end



clear T_reg V_reg RHfp_reg RHhumi_reg RHfp_i_reg RHhumi_i_reg LWCfp_reg LWChumi_reg Displacement_reg qv_humi_reg
for i=1:length(x)
    T_reg(:,i)=interp1(YY(:,i),TT(:,i),z);
    V_reg(:,i)=interp1(YY(:,i),VV(:,i),z);
    RHfp_reg(:,i)=interp1(YY(:,i),RH_fp_ALL(:,i),z);
    RHhumi_reg(:,i)=interp1(YY(:,i),RH_humi_ALL(:,i),z);  
    RHfp_i_reg(:,i)=interp1(YY(:,i),RH_fp_i_ALL(:,i),z);
    RHhumi_i_reg(:,i)=interp1(YY(:,i),RH_humi_i_ALL(:,i),z); 
    LWCfp_reg(:,i)=interp1(YY(:,i),LWC_fp_ALL(:,i),z);
    LWChumi_reg(:,i)=interp1(YY(:,i),LWC_humi_ALL(:,i),z);   
    Displacement_reg(:,i)=interp1(YY(:,i),YY(:,i)-YY(:,1),z);
    pot_reg(:,i)=interp1(YY(:,i),pot_wave(:,i),z);
    qv_humi_reg(:,i)=interp1(YY(:,i),qv_wave_humi(:,i),z);
end

%P_reg = interp1(h,P_wave,z);
P_reg = repmat(P_wave,[1 size(T_reg,2)]);
%pot_reg = (T_reg+273.15).* (1e5./P_reg).^0.286;
% 
% ZZ=100 + 0.35*y_upper; %create a shallow sin path through for the aircraft
% T_path=interp2(x,z,T_reg,x,ZZ); %interpolate from the 2D field along the path

airplane_dis = 75*sin(x*xf); %give the aircraft the same wavelenght with an amplitude of 150 m as observed
%y=0 displacement of plane
plane_z0 = z0+yf*1000; %through the crest of the big wave
plane_z0 = 3150; %or choose a heigt

%moved all these down by 100 cf. original
ZZ=+100 + airplane_dis + plane_z0; %create a shallow sin path through for the aircraft
T_path=interp2(x,z,T_reg,x,ZZ); %adding z0+plane_z0 to centre them at the crest of the big wave
V_path=interp2(x,z,V_reg,x,ZZ);
LWChumi_path=interp2(x,z,LWChumi_reg,x,ZZ);

ZZ2=+50 + airplane_dis + plane_z0; %create a shallow sin path through for the aircraft
T_path2=interp2(x,z,T_reg,x,ZZ2);
V_path2=interp2(x,z,V_reg,x,ZZ2);
LWChumi_path2=interp2(x,z,LWChumi_reg,x,ZZ2);

%ZZ2=100 + 0.25*y_upper; %create a shallow sin path through for the aircraft
%T_path2=interp2(x,z,T_reg,x,ZZ2); %interpolate from the 2D field along the path

%have chosen ZZ3 to be of the same amplitude as T_upper3 as this
%had an amplitude of 1.5 degrees, similar to what the temperature
%oscillation would be for air that moved at the same amplitude as that
%which the aircraft moved (based on the pressure of the aircraft)
ZZ3=-50 + airplane_dis+plane_z0; %create a shallow sin path through for the aircraft
T_path3=interp2(x,z,T_reg,x,ZZ3); %interpolate from the 2D field along the path
V_path3=interp2(x,z,V_reg,x,ZZ3);
LWChumi_path3=interp2(x,z,LWChumi_reg,x,ZZ3);

%same, but higher up
ZZ4= 0 + airplane_dis+plane_z0; %create a shallow sin path through for the aircraft
T_path4=interp2(x,z,T_reg,x,ZZ4); %interpolate from the 2D field along the path
V_path4=interp2(x,z,V_reg,x,ZZ4);
LWChumi_path4=interp2(x,z,LWChumi_reg,x,ZZ4);

append='_idamp=1';
%savedir='C:\Documents and Settings\dan\My Documents\MATLAB\Antarctica\detailed_look_at_flights\lee_wave_case\';


X=repmat(x,[size(RHhumi_reg,1) 1]);
Z=repmat(z,[size(RHhumi_reg,2) 1])';

clear T_max qvsat_max qv_max equiv_max_fp

[ymax imax]=max(y);
P_max=interp1(h,P_wave,z);
for i=1:length(imax)
    T_max(i)=273.15+T_reg(i,imax(i));
    qv_max(i)=SatVapPress(T_max(i),'goff','liq',P_max(i),1)/f; %assume saturation
    equiv_max_fp(i) = equivalent_potemp(T_max(i),P_max(i),qv_max(i));
end

%P_reg = repmat(P_max,[size(T_reg,2) 1])';
for i=1:size(T_reg,2)
    qsat_reg = SatVapPress(273.15+T_reg(:,i)','goff','liq',P_max,1)/f; %assume saturation
    equiv_reg_humi(:,i) = equivalent_potemp(273.15+T_reg(:,i)',P_max,qsat_reg);
end

% equiv_reg_humi(icut_off:end,:)=NaN;

cutout_low = 3075;
cutout_high = 3400; %cuts out the values between these heights.
icut_out = find(z<cutout_low);
icut_out = [icut_out find(z>cutout_high)];
%equiv_max_fp(icut_off:end)=NaN;

%want to find the max LWC for each equiv potemp. Will do it in bins of equiv potemp
%as the values are not regular
bins=[minALL(equiv_reg_humi):0.2:maxALL(equiv_reg_humi)];
E=equiv_reg_humi;
L=LWChumi_reg;
E(icut_out,:)=''; %cut off values with the height above 3400 m
L(icut_out,:)='';
[meanvals,equiv_max_humi2,max_inds]=bin_data(E(:),L(:),bins); %finds the max and mean of the bins



if iplot==1

tit='Temperature (^{o}C) field of gravity waves with aircraft paths marked';
figure('name',tit);
pcolor(x,z,T_reg); shading interp; colorbar; hold on
plot(x,ZZ,'k');
plot(x,ZZ2,'k');
plot(x,ZZ3,'k--');
plot(x,ZZ4,'r--');
%plot(x,y_lower+z0);
%plot(x,y+z0);
%plot(x,y_upper+z0);
%plot(x,y_upper2+z0);
%plot(x,y_upper3+z0);
title(tit);
xlabel('X (km)');
ylabel('Z (m)');
if isave==1
    saveas(gcf,[filedir tit append '.fig'],'fig');
    print(gcf,[filedir tit append '.emf'],'-dmeta');    
    close    
end

tit='Vertical wind speed (m s^{-1}) field of gravity waves with aircraft paths marked';
figure('name',tit);
pcolor(x,z,V_reg); shading interp; colorbar; hold on
plot(x,ZZ,'k');
plot(x,ZZ2,'k');
plot(x,ZZ3,'k--');
plot(x,ZZ4,'r--');
% plot(x,y_lower+z0);
% plot(x,y+z0);
% plot(x,y_upper+z0);
% plot(x,y_upper2+z0);
% plot(x,y_upper3+z0);
title(tit);
xlabel('X (km)');
ylabel('Z (m)');
if isave==1
    saveas(gcf,[filedir tit append '.fig'],'fig');
    print(gcf,[filedir tit append '.emf'],'-dmeta');    
    close    
end

tit='Hygrometer instrument RH (%) field of gravity waves with aircraft paths marked';
figure('name',tit);
pcolor(x,z,RHfp_reg); shading interp; colorbar; hold on
plot(x,ZZ,'k');
plot(x,ZZ2,'k');
plot(x,ZZ3,'k--');
plot(x,ZZ4,'r--');
% plot(x,y_lower+z0);
% plot(x,y+z0);
% plot(x,y_upper+z0);
% plot(x,y_upper2+z0);
% plot(x,y_upper3+z0);
title(tit);
xlabel('X (km)');
ylabel('Z (m)');

i=find(abs(RHfp_reg-100)<0.1);
plot(X(i),Z(i),'ko','markerfacecolor','k','markersize',2);

if isave==1
    saveas(gcf,[filedir tit append '.fig'],'fig');
    print(gcf,[filedir tit append '.emf'],'-dmeta');    
    close    
end

tit='Hygrometer instrument RHi (%) field of gravity waves with aircraft paths marked';
figure('name',tit);
pcolor(x,z,RHfp_i_reg); shading interp; colorbar; hold on
plot(x,ZZ,'k');
plot(x,ZZ2,'k');
plot(x,ZZ3,'k--');
plot(x,ZZ4,'r--');
% plot(x,y_lower+z0);
% plot(x,y+z0);
% plot(x,y_upper+z0);
% plot(x,y_upper2+z0);
% plot(x,y_upper3+z0);
title(tit);
xlabel('X (km)');
ylabel('Z (m)');

i=find(abs(RHfp_i_reg-100)<0.1);
plot(X(i),Z(i),'ko','markerfacecolor','k','markersize',2);

if isave==1
    saveas(gcf,[filedir tit append '.fig'],'fig');
    print(gcf,[filedir tit append '.emf'],'-dmeta');    
    close    
end

tit='Humicap RH (%) field of gravity waves with aircraft paths marked';
figure('name',tit);
pcolor(x,z,RHhumi_reg); shading interp; colorbar; hold on
plot(x,ZZ,'k');
plot(x,ZZ2,'k');
plot(x,ZZ3,'k--');
plot(x,ZZ4,'r--');
% plot(x,y_lower+z0);
% plot(x,y+z0);
% plot(x,y_upper+z0);
% plot(x,y_upper2+z0);
% plot(x,y_upper3+z0);
title(tit);
xlabel('X (km)');
ylabel('Z (m)');

%contour(x,z,RHhumi_reg,[100 100],'color','g','linewidth',2);
i=find(abs(RHhumi_reg-100)<0.1);
plot(X(i),Z(i),'ko','markerfacecolor','k','markersize',2);
if isave==1
    saveas(gcf,[filedir tit append '.fig'],'fig');
    print(gcf,[filedir tit append '.emf'],'-dmeta');    
    close    
end

tit='Humicap RHi (%) field of gravity waves with aircraft paths marked';
figure('name',tit);
pcolor(x,z,RHhumi_i_reg); shading interp; colorbar; hold on
plot(x,ZZ,'k');
plot(x,ZZ2,'k');
plot(x,ZZ3,'k--');
plot(x,ZZ4,'r--');
% plot(x,y_lower+z0);
% plot(x,y+z0);
% plot(x,y_upper+z0);
% plot(x,y_upper2+z0);
% plot(x,y_upper3+z0);
title(tit);
xlabel('X (km)');
ylabel('Z (m)');
i=find(abs(RHhumi_i_reg-100)<0.1);
plot(X(i),Z(i),'ko','markerfacecolor','k','markersize',2);
if isave==1
    saveas(gcf,[filedir tit append '.fig'],'fig');
    print(gcf,[filedir tit append '.emf'],'-dmeta');    
    close    
end

tit='Hygrometer LWC (g m^{-3}) field of gravity waves with aircraft paths marked';
figure('name',tit);
pcolor(x,z,LWCfp_reg); shading interp; colorbar; hold on
plot(x,ZZ,'k');
plot(x,ZZ2,'k');
plot(x,ZZ3,'k--');
plot(x,ZZ4,'r--');
% plot(x,y_lower+z0);
% plot(x,y+z0);
% plot(x,y_upper+z0);
% plot(x,y_upper2+z0);
% plot(x,y_upper3+z0);
title(tit);
xlabel('X (km)');
ylabel('Z (m)');
%i=find(abs(RHhumi_i_reg-100)<0.1);
%plot(X(i),Z(i),'ko','markerfacecolor','k','markersize',2);
if isave==1
    saveas(gcf,[filedir tit append '.fig'],'fig');
    print(gcf,[filedir tit append '.emf'],'-dmeta');    
    close    
end

tit='Humicap LWC (g m^{-3}) field of gravity waves with aircraft paths marked';
figure('name',tit);
pcolor(x,z,LWChumi_reg); shading interp; colorbar; hold on
plot(x,ZZ,'k');
plot(x,ZZ2,'k');
plot(x,ZZ3,'k--');
plot(x,ZZ4,'r--');
% plot(x,y_lower+z0);
% plot(x,y+z0);
% plot(x,y_upper+z0);
% plot(x,y_upper2+z0);
% plot(x,y_upper3+z0);
title(tit);
xlabel('X (km)');
ylabel('Z (m)');
%i=find(abs(RHhumi_i_reg-100)<0.1);
%plot(X(i),Z(i),'ko','markerfacecolor','k','markersize',2);
if isave==1
    saveas(gcf,[filedir tit append '.fig'],'fig');
    print(gcf,[filedir tit append '.emf'],'-dmeta');
    close    
end


tit='Temperature change along aircraft paths';
figure('name',tit);
plot(x,T_path,'b'); set(gca,'ydir','reverse'); hold on
plot(x,T_path2,'r');
plot(x,T_path3,'r--');
plot(x,T_path4,'g--');
ylabel('T ^{o}C');
xlabel('X (km)');
title(tit);
if isave==1
    saveas(gcf,[filedir tit append '.fig'],'fig');
    print(gcf,[filedir tit append '.emf'],'-dmeta');    
    close    
end

tit='Vertical velocity change along aircraft paths';
figure('name',tit);
plot(x,V_path,'b'); hold on
plot(x,V_path2,'r');
plot(x,V_path3,'r--');
plot(x,V_path4,'g--');
ylabel('W (m s^{-1})');
xlabel('X (km)');
title(tit);
if isave==1
    saveas(gcf,[filedir tit append '.fig'],'fig');
    print(gcf,[filedir tit append '.emf'],'-dmeta'); 
    close    
end

tit='LWC along aircraft paths';
figure('name',tit);
plot(x,LWChumi_path,'b'); hold on
plot(x,LWChumi_path2,'r');
plot(x,LWChumi_path3,'r--');
plot(x,LWChumi_path4,'g--');
ylabel('LWC (g m^{-3})');
xlabel('X (km)');
title(tit);
if isave==1
    saveas(gcf,[filedir tit append '.fig'],'fig');
    print(gcf,[filedir tit append '.emf'],'-dmeta');    
    close    
end


tit='Hygrometer LWC profile at max displacement';
tit='Max hygrometer LWC profile';
figure('name',tit);
%plot(LWCfp_reg(:,imax),z,'r'); hold on
[max_LWC_fp,imax]=max(LWCfp_reg,[],2);
plot(max_LWC_fp,z,'r'); hold on


ylabel('Z (m)');
xlabel('LWC g m^{-3}');
title(tit);
set(gca,'ylim',[2500 4000]);
if isave==1
    saveas(gcf,[filedir tit append '.fig'],'fig');
    print(gcf,[filedir tit append '.emf'],'-dmeta'); 
    close    
end

%tit='Humicap LWC profile at max displacement';
tit='Max humicap LWC profile';
figure('name',tit);

[max_LWC_humi,imax]=max(LWChumi_reg,[],2);
plot(max_LWC_humi,z,'r'); hold on

clear T_max qvsat_max qv_max equiv_max_humi equiv_reg_humi
P_max=interp1(h,P_wave,z);
for i=1:length(imax)
    T_max(i)=273.15+T_reg(i,imax(i));
    qv_max(i)=SatVapPress(T_max(i),'goff','liq',P_max(i),1)/f; %assume saturation
    equiv_max_humi(i) = equivalent_potemp(T_max(i),P_max(i),qv_max(i));
end






% 
% [max_eq_humi,imax]=max(equiv_reg_humi,[],2);
% %clear T_max qvsat_max qv_max equiv_max_humi equiv_reg_humi
% for i=1:i
%     T_max(i)=273.15+T_reg(i,imax(i));
%     qv_max(i)=SatVapPress(T_max(i),'goff','liq',P_max(i),1)/f; %assume saturation
%     equiv_max_humi(i) = equivalent_potemp(T_max(i),P_max(i),qv_max(i));
% end


ylabel('Z (m)');
xlabel('LWC g m^{-3}');
title(tit);
set(gca,'ylim',[2500 3500]);
if isave==1
    saveas(gcf,[filedir tit append '.fig'],'fig');
    print(gcf,[filedir tit append '.emf'],'-dmeta');  
    close    
end


tit='Hotwire scatter plot overlay';
ichoose_scatter_vert_coord=1;
vertical_coord='height'; 
scatter_plot
plot(max(LWChumi_reg,[],2),z/1000,'r','linewidth',2); hold on
plot(max(LWCfp_reg,[],2),z/1000,'k','linewidth',2); hold on
if isave==1
    saveas(gcf,[filedir tit append '.fig'],'fig');
    print(gcf,[filedir tit append '.emf'],'-dmeta');    
    close    
end


tit='Hotwire vs equiv potemp scatter plot overlay';
ichoose_scatter_vert_coord=1;
vertical_coord='equiv potemp sat'; 
scatter_plot
%plot(max_LWC_humi,equiv_max_fp,'r','linewidth',2); hold on
plot(equiv_max_humi2,0.5*(bins(1:end-1)+bins(2:end)),'r','linewidth',2);
%equiv_max_humi2 is LWC, bins are equiv potemp

if isave==1
    saveas(gcf,[filedir tit append '.fig'],'fig');
    print(gcf,[filedir tit append '.emf'],'-dmeta');  
    close    
end


if iprofs==1
    tit='Water Vapour profile';
    figure    
    plot(qfp_sound,Zsound,'linewidth',2); hold on
    plot(qhumi_sound,Zsound,'r','linewidth',2); hold on
    plot(vap_factor*qv_prof_fp,h,'b--','linewidth',2);
    plot(vap_factor*qv_prof_humi,h,'r--','linewidth',2);
    grid
    set(gca,'ylim',[2000 3500]);
    xlabel('Qv (kg kg^{-1})');
    ylabel('Height (m)');
    clear leg_str
    leg_str{1}='Sounding FP';
    leg_str{2}='Sounding humi';
    leg_str{3}='Idealized FP';
    leg_str{4}='Idealized humi';  
    legend(leg_str,'Location','NorthWest');

    
    if isave==1
        saveas(gcf,[filedir tit append '.fig'],'fig');
        print(gcf,[filedir tit append '.emf'],'-dmeta');    
        close        
    end

    
    
    tit='Temperature profile';
    figure    
    plot(Tsound,Zsound,'linewidth',2); hold on
    plot(T_h-273.15,h,'b--','linewidth',2);
    grid
    set(gca,'ylim',[2000 3500]);
    xlabel('Temperature (^{o}C)');
    ylabel('Height (m)');
    clear leg_str
    leg_str{1}='Sounding';
    leg_str{2}='Idealized';
    legend(leg_str,'Location','NorthEast');

    
    if isave==1
        saveas(gcf,[filedir tit append '.fig'],'fig');
        print(gcf,[filedir tit append '.emf'],'-dmeta');   
        close
    end


    tit='Equivalent Potential Temperature profile';
    figure    
    plot(equiv_sound_fp,Zsound,'linewidth',2); hold on
    plot(equiv_sound_humi,Zsound,'r','linewidth',2); hold on
    plot(eq_h_fp,h,'b--','linewidth',2);
    plot(eq_h_humi,h,'r--','linewidth',2);
    grid
    set(gca,'ylim',[2000 3500]);
    xlabel('Equiv. Potemp (K)');
    ylabel('Height (m)');
    clear leg_str
    leg_str{1}='Sounding FP';
    leg_str{2}='Sounding humi';
    leg_str{3}='Idealized FP';
    leg_str{4}='Idealized humi';  
    legend(leg_str,'Location','NorthWest');

    
    if isave==1
        saveas(gcf,[filedir tit append '.fig'],'fig');
        print(gcf,[filedir tit append '.emf'],'-dmeta');    
        close
    end
    



    
    
end


end


if length(z0damp)>1
    z0damp=0;
end

fprintf(1,'\nidamp=%.1f, vap_factor=%.2f, qvfp_mean=%2.1e, z0damp=%.0f m, dz_damp=%.0f m, iuse_equiv_constant_lapse_rate=%d, iuse_constant_temp_lapse_rate=%d, add_constant_temp=%1.2f, add_constant_eq=%1.2f, yf=%3.0f m, iuse_constant_RH_above_sounding=%d, RH=%1.2f\n'...
    ,idamp,vap_factor,qvfp_mean,z0damp,dz_damp,iuse_equiv_constant_lapse_rate,iuse_constant_temp_lapse_rate,...
    add_constant_temp,add_constant_eq,yf*1000,iuse_constant_RH_above_sounding,RH_constant);






