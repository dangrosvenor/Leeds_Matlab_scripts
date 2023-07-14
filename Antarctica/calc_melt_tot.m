%icorrect_alb_emiss=0;
icorrect_alb_emiss=1; %Flag to fix the albedo and emissivity to the aircraft observed values
%alb = 0.78, emiss=1.0

melt_tot=zeros(size(lat2d.var));
melt_tot_total=zeros(size(lat2d.var));
sw_tot=zeros(size(lat2d.var));
lw_tot=zeros(size(lat2d.var));
lwdown_tot=zeros(size(lat2d.var));
sh_tot=zeros(size(lat2d.var));
lh_tot=zeros(size(lat2d.var));
grd_tot=zeros(size(lat2d.var));
temp_tot=zeros(size(lat2d.var));
sp10_tot=zeros(size(lat2d.var));
cond_tot=zeros(size(lat2d.var));
n_tot=zeros(size(lat2d.var));
n_tot2=zeros(size(lat2d.var));
rh_tot=zeros(size(lat2d.var));
rh_iz_tot=zeros(size(lat2d.var));
max_melt=zeros(size(lat2d.var));
max_melt_time=zeros(size(lat2d.var));
tsk_min=NaN*ones(size(lat2d.var));


ih_wrf=1; %height level for RH averages

same_times='yes';     %flag to say whether to use the same times for each location regardless of whether the skin temp is zero or not 
                      %(i.e. regardless of if is melting)
%same_times='no';     %or use any times but only when melting is occurring (skin temp=0)

only_one_day=0; %flag to say that we only want to calculate data for a particular 
only_one_day=6; %day. Set to 0 for all days or to N for date N (e.g. set to
%6 for 6th Jan (ignoring the month here).

av_type='total'; %add up all the contribution values for a 3 day total mean
%av_type='mean';  %take average of the contributions only for the times when melting was occurring for each point

sample_times = [12 15 18 21]; %times to use for the samples (if same_times='yes')
  %is always melting over the ice shelves at these times - but outside of
  %these times the skin temp does drop below zero - so need to put in
  %check for TSK>=0 only
%sample_times = [0 3 6 9]; %times to use for the samples (if same_times='yes')
%sample_times = [12];

f=1e6*28.97/18; %conversion between MR and ppmv - use 18 for water vapour and 48 for ozone

clear hours_num
switch same_times
              case 'yes'
                  hours_str=Times(:,12:13); %text of the hours in wrfout
                  days_str=Times(:,9:10); %text of the hours in wrfout                  
                  for i=1:size(hours_str,1)
                      hours_num(i) = str2num(hours_str(i,:)); %convert the two rows of characters to a 1D number vector
                      days_num(i) = str2num(days_str(i,:)); %convert the two rows of characters to a 1D number vector                      
                  end

                  for itime=1:length(sample_times)
                      i=find(hours_num==sample_times(itime)); %find the indices of the required times in Times
                          if only_one_day~=0
                              i=find(hours_num==sample_times(itime) & days_num==only_one_day)
                          end
                      if itime==1
                          times = i;
                      else
                          times = [times i];
                      end
                  end
    otherwise
        times=2:25;
%        times=9;
%        times=[1     9    17    25     2    10    18     3    11    19     4    12    20]; %indices for times other than 12,15,18&21
%         times=[5    13    21     6    14    22     7    15    23     8    16    24]; %12,15,18,21 UTC for 3 days
end

time_step = diff(days_num*24+hours_num);
if max(time_step)~=min(time_step)
    disp('*** Problem with the time_step !!! ***');
else
    time_step=time_step(1);
end


for time=times
        time
        LW=nc{'GLW'}(time,:);  %downwelling LW at surface
        SW=nc{'SWDOWN'}(time,:); %downwelling SW at surface
        SH=-nc{'HFX'}(time,:);
        LH=-nc{'LH'}(time,:); %negative as WRF convention is that these are fluxes into air
        ALBEDO=nc{'ALBEDO'}(time,:);
        EMISS=nc{'EMISS'}(time,:);
        GRDFLX=nc{'GRDFLX'}(time,:); %sign convention is that this is flux into the surface layer
        TSK=nc{'TSK'}(time,:);
        T2=nc{'T2'}(time,:);
        u10=nc{'U10'}(time,:);
        v10=nc{'V10'}(time,:);
        sp10=sqrt(u10.^2+v10.^2);
        q2=nc{'Q2'}(time,:);
        psfc=nc{'PSFC'}(time,:);
        q_iz=nc{'QVAPOR'}(time,ih_wrf,:);
        P_iz=nc{'P'}(time,ih_wrf,:) + nc{'PB'}(time,ih_wrf,:);
        T_iz = (nc{'T'}(time,ih_wrf,:) + 300) ./ ( (1e5./P_iz).^0.286 );
       
        
        if icorrect_alb_emiss==1
            ALBEDO = 0.78;  %WRF value is 0.7
        end
        
        %sign convention is postive means energy going into ground
        SW_UP=ALBEDO.*SW; %think albedo means this for WRF
        SW_NET=SW-SW_UP;

        %EMISS=0.95; %Polar WRF value should be set to 0.98 - however, the
        %output says that it is 0.98 for the first time step, but then it
        %changes to 0.95 for all subsequent ones. Is this just a bug in the
        %output, or for the calculations themselves?
        
%        LW_UP=5.67e-8.*TSK.^4; %Boltzmann W/m2 using skin temperature
%        LW_NET=EMISS.*(LW-LW_UP); %Should be LW_NEW = LW - EMISS.*LW_UP; ??Since emissivity just affects the surface upwelling?
          %and not the downwelling (LW)?
          
          if icorrect_alb_emiss==1
              EMISS=1.0; %aircraft derived value was almost one
              fprintf(1,'\nWARNING - fixing model emissivity and albedo to aircraft values');
          end
          
          %Or better as
         LW_UP=EMISS*5.67e-8.*TSK.^4; %Boltzmann W/m2 using skin temperature 
         LW_NET = LW-LW_UP; %LW is the downwelling LW at surface, so don't use emissivity for this
          
          
        
        landmask=nc{'LANDMASK'}(1,:);
        seaice=nc{'SEAICE'}(1,:); %NOTE seaice is zero everywhere for NCEP as it wasn't put in properly on this run!
                                  %think have a run with seaice included

        rh = q2./(satvappress(T2,'goff','liq',psfc,1)/f); %using the surface pressure to calc RH
        rh(rh>1)=1;
        
        rh_iz = q_iz./(satvappress(T_iz,'goff','liq',P_iz,1)/f); %using the surface pressure to calc RH
        rh_iz(rh_iz>1)=1;

          switch same_times
              case 'yes'
                  i0=find( abs(landmask-1)<0.001 & abs(seaice-0)<0.000001); %find all points that are land (but not sea ice)
                  tsk_temp =NaN*ones(size(lat2d.var));
                  tsk_temp(i0) = TSK(i0)-273.15;
                  tsk_min = min(tsk_min,tsk_temp);
                  i0=find( abs(landmask-1)<0.001 & abs(seaice-0)<0.000001 & TSK>273.15-0.001); %find all points that are land (but not sea ice)                  
                  i02=find( abs(landmask-1)<0.001 & abs(seaice-0)<0.000001 & abs(TSK-273.15)>0.001 ); %land, but not seaice and surf temp not zero deg C
              case 'no'
                  i0=find( abs(TSK-273.15)<0.001 & abs(landmask-1)<0.001 & abs(seaice-0)<0.000001); %find all points where skin temp was zero deg C and that are land
                  i02=find( abs(TSK-273.15)>0.001 & abs(landmask-1)<0.001 & abs(seaice-0)<0.000001); %find all points where skin temp wasn't zero deg C and that are land
          end
                                  

        
        
        i_condensate=0;
        recalc=1;
        if i_condensate==1
            if recalc==1
                ih_inds=1:30;

                M = 28.97*1.67E-27;
                k = 1.38E-23;
                P = WRFUserARW(nc,'p',time); P.var=P.var*100;
                T = WRFUserARW(nc,'tc',time); T.var=T.var+273.15;
                rho=P.var(ih_inds,:,:).*M./k./T.var(ih_inds,:,:);
                Z = squeeze(WRFUserARW(nc,'Z',time));
                dz_grid = Z.var(ih_inds+1,:,:)-Z.var(ih_inds,:,:);

%                cond_dat = squeeze (sum( (nc{'QCLOUD'}(time,ih_inds,:,:) + nc{'QRAIN'}(time,ih_inds,:,:)...
%                    + nc{'QICE'}(time,ih_inds,:,:) + nc{'QSNOW'}(time,ih_inds,:,:)...
%                    + nc{'QGRAUP'}(time,ih_inds,:,:)) .* rho*1000*dx_grid*1000*dy_grid.*dz_grid ,1) );
                
%not multiplying by dx and dy to keep as kg/m2 so that is normalised by area - otherwise will need to know the grid
%cell size for the value to mean anything.
                cond_dat = squeeze (sum( (nc{'QCLOUD'}(time,ih_inds,:,:) + nc{'QRAIN'}(time,ih_inds,:,:)...
                    + nc{'QICE'}(time,ih_inds,:,:) + nc{'QSNOW'}(time,ih_inds,:,:)...
                    + nc{'QGRAUP'}(time,ih_inds,:,:)) .* rho.*dz_grid ,1) );
                
            else
                fprintf(1,'\n****** WARNING - not recalculating melt energy values *******');
            end                                    

        end
        

        switch av_type
            case 'total'
                
                %        i0=1:prod(size(melt_tot));  %for all points
                %fluxes are in W/m2
                melt_tot(i0) = melt_tot(i0) + ( SW_NET(i0) + LW_NET(i0) + LH(i0) + SH(i0) + GRDFLX(i0) );

                %multiply by 3 (the timestep of the runs)
                %and divide by 72 to give the average over the 3 day period with data every 3 hours
                %convert to per day and divide by the latent heat of fusion to give in mm/day
        
                sw_tot(i0) = sw_tot(i0)+ SW_NET(i0);  %for the mean across the whole 72 hours (so lower number of times when T=0 will reduce the mean)
                lw_tot(i0) = lw_tot(i0)+ LW_NET(i0);  %mean then applies for whole 3 days so total melting is this multiplied by 3
                lwdown_tot = lwdown_tot + LW;
                sh_tot(i0) = sh_tot(i0) + SH(i0);
                lh_tot(i0) = lh_tot(i0) + LH(i0);
                grd_tot(i0) = grd_tot(i0) + GRDFLX(i0);

            case 'mean'
                melt_tot_total(i0) = melt_tot_total(i0) + ( SW_NET(i0) + LW_NET(i0) + LH(i0) + SH(i0) + GRDFLX(i0) );
                melt_tot(i0) = melt_tot(i0) + SW_NET(i0) + LW_NET(i0) + LH(i0) + SH(i0) + GRDFLX(i0) ;
                sw_tot(i0) = sw_tot(i0)+ SW_NET(i0);  %for the mean across the whole 72 hours (so lower number of times when T=0 will reduce the mean)
                lw_tot(i0) = lw_tot(i0)+ LW_NET(i0);  %mean then applies for whole 3 days so total melting is this multiplied by 3
                lwdown_tot = lwdown_tot + LW;
                sh_tot(i0) = sh_tot(i0) + SH(i0);
                lh_tot(i0) = lh_tot(i0) + LH(i0);
                grd_tot(i0) = grd_tot(i0) + GRDFLX(i0);


        end
        
        max_melt_t = zeros(size(SW_NET));
        max_melt_t(i0) = SW_NET(i0) + LW_NET(i0) + LH(i0) + SH(i0) + GRDFLX(i0);
        
        imax_t = find(max_melt_t>max_melt);
        max_melt(imax_t)=max_melt_t(imax_t);
        max_melt_time(imax_t)=time;


        sp10_tot(i0)= sp10_tot(i0) + sp10(i0);    %for the mean within the samples we have (when T=0)
        temp_tot(i0)= temp_tot(i0) + T2(i0);
        rh_tot(i0) = rh_tot(i0) + rh(i0);
        rh_iz_tot(i0) = rh_iz_tot(i0) + rh_iz(i0);
        
        if i_condensate==1
            cond_tot(i0)= cond_tot(i0) + cond_dat(i0);
        end
        n_tot(i0) = n_tot(i0) + 1;
        n_tot2(i02) = n_tot2(i02) + 1; 


end


%melt_tot etc are in W/m2 (J/s/m2) at the moment
%latent heat of fusion = 3.34e5 J/kg. time_step*3600 the number of seconds in one
%timestep
%Divide by latent heat to convert to mass of water melted per unit area in kg/m2.
%1 mm of water per m2 has a mass of 1 kg. So 1 kg/m2 is equivalent to 1 mm
%of depth. (density of water=1000 kg/m3, so 1 mm x 1m2 is a volume 1e-3*1
%m3 and so weighs 1 kg.

melt_tot=melt_tot*time_step*3600 / 3.34e5;  %convert to total melt contribution in mm (total over 72 hours)
sw_tot=sw_tot*time_step*3600 / 3.34e5;
lw_tot=lw_tot*time_step*3600 / 3.34e5;
lwdown_tot=lwdown_tot*time_step*3600 / 3.34e5;
sh_tot=sh_tot*time_step*3600 / 3.34e5;
lh_tot=lh_tot*time_step*3600 / 3.34e5;
grd_tot=grd_tot*time_step*3600 / 3.34e5;

temp_tot=temp_tot./n_tot;  %divide by number of samples to get the mean within the samples
sp10_tot=sp10_tot./n_tot;  
cond_tot=cond_tot./n_tot;  
rh_tot=rh_tot./n_tot;
rh_iz_tot=rh_iz_tot./n_tot;

switch av_type
    case 'mean'
        melt_tot_total = melt_tot_total*time_step*3600 / 3.34e5; %convert to the total melting that occurs in the 3 days (mm)
        melt_tot=melt_tot./n_tot/3;  %convert to the mean over just the melting times
        sw_tot=sw_tot./n_tot/3;  %this then gives the average melting rate for the times when melting occurred in mm/day
        lw_tot=lw_tot./n_tot/3;
        lwdown_tot=lwdown_tot./n_tot/3;
        sh_tot=sh_tot./n_tot/3;
        lh_tot=lh_tot./n_tot/3;
        grd_tot=grd_tot./n_tot/3;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    calculate the means along each latitude ***************
heat_fluxes_mean_along_latitude %calculate the mean values along each latitude


