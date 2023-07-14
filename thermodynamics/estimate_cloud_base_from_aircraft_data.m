%estimate cloud base from the aircraft data

%first of all calculated mean profiles from the aircraft data by binning in z

dz=10; %dz in m - probably don't want to make this too small as would
%then get large fluctuations with height

z0=400; %estimate a height for the starting qv (perhaps should use an average?)
  
%OR could use method where we set the air parcel start conditions and then calculate cloud
%base pressure from that
q_parcel = 2.9e-3; %mixing ratio kg/kg
p_parcel = 990; %set the starting pressure
th_parcel = 263.79; %set the starting potemp
    
    
     times=[13.25 14];   
     times=[];
     times=[16.0 time_flt(end)]; 
        
clear indsALL
        if length(times)==0
            indsALL(1).i=1:length(dat_flt(:,1));
        else
            for i_inds=1:size(times,1)
                [time_0,time_1]=findheight_nearest(dat_flt(:,1)/1000/3600,times(i_inds,1),times(i_inds,2));
                indsALL(i_inds).i=time_0:time_1;
            end
        end

        inds=indsALL(1).i;
        
        
    Z=dat_flt(inds,col_alt);
    T=dat_flt(inds,col_temp);
    P=dat_flt(inds,col_press); 
    qv=qv_flt100_humi(inds);
%    qv=qv_flt100_fp(inds);    
        
    
    clear tgrid pgrid qvgrid
    zgrid=[0:dz:max(Z)];
    for iz=1:length(zgrid)-1
        ii=find( Z>zgrid(iz) & Z<zgrid(iz+1) );
        tgrid(iz)=mean(T(ii));
        pgrid(iz)=mean(P(ii));
        qvgrid(iz)=mean(qv(ii));
        if length(ii)==0
            tgrid(iz)=NaN;
        end
    end

    zgrid=zgrid(2:end);
    
    inan=isnan(tgrid);
    inan2=find(inan==1);
    tgrid(inan2)=[];
    pgrid(inan2)=[];
    qvgrid(inan2)=[];
    zgrid(inan2)=[];
    
    pot_grid=(273.15+tgrid).*(1000./pgrid).^0.286;
    
    

    th=interp1(zgrid,pot_grid,z0);
    p=interp1(zgrid,pgrid,z0);
    q=interp1(zgrid,qvgrid,z0);
    t=interp1(zgrid,tgrid,z0);    
    
    Tad=th./(1000./pgrid).^0.286; %calculate the temperatures during a constant dry adiabatic ascent
    
    t_dew=Tdew(q,p*100)-273.15; %the dew point temperature for our q
    %what about the pressure variation with height?
    
    %find the location in Tad that is closest to t_dew

    
    [Tad_uni iuni]=unique(Tad);
    z_cb = interp1(Tad_uni,zgrid(iuni),t_dew+273.15)
    
    %% alternative method where the parcel base conditions are set (not using the aircraft data
    %% except to find the height for a given cloud base pressure
    pgrid_parcel = [1080:-1:100];
    Tad_parcel=th_parcel./(1000./pgrid_parcel).^0.286; %calculate the temperatures during a constant dry adiabatic ascent    
    t_dew_parcel=Tdew(q_parcel,p_parcel*100)-273.15; %the dew point temperature for our q
    [minT_parcel imin_parcel]=min(abs(Tad_parcel-273.15-t_dew_parcel));
    cb_pressure_parcel = pgrid_parcel(imin_parcel)
    if cb_pressure_parcel>pgrid(1)
        disp('Cloud base is below the max pressure in the aircraft data');
    else
        cb_height_parcel =interp1(pgrid,zgrid,cb_pressure_parcel)
    end
    
    % *** alternative method ***
    %here will iterate over height re-calculating dew point each time
    min_diff=1e9;
    for it=1:length(Tad)
        t_dew2=Tdew(q,pgrid(it)*100);
        dif=abs(t_dew2-Tad(it));
        if dif<min_diff
            it_cb=it;
            min_diff=dif;
        end
        
    end
    
    z_cb2=zgrid(it_cb)
    
    
    
    
    
    