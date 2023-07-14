% 
% *** clearing here ***
clear modis_var; istring=1;
modis_var{istring}='Soil_Moisture'; istring=istring+1; %

MLAT = nc{'lat'}(:);
MLON = nc{'lon'}(:);


for iread_modis=1:length(modis_var)

    MODIS_varname = modis_var{iread_modis};

    %get the index number from the name
%    nvar = hdfsd('nametoindex',SD_id,MODIS_varname)+1; %(nvar=316);
    
    %read the data
%    modis_data_read = double(hdfread(INFO.Vgroup(1).Vgroup(1).SDS(nvar))); %this retrieves the data   
    var_nc = nc{MODIS_varname};
    modis_data_read = var_nc(:);
    
    
    
    if length(size(modis_data_read))==4   %position of the units, scalef_factor and offset change
        iunits=6;                    %depending on the dimension of the array
        Y_hbins = double(INFO.Vgroup(1).Vgroup(1).SDS(nvar).Attributes(3).Value);
        X_hbins = double(INFO.Vgroup(1).Vgroup(1).SDS(nvar).Attributes(4).Value);     
        eval([MODIS_varname '.Xbins=X_hbins;']);
        eval([MODIS_varname '.Ybins=Y_hbins;']);
        
        ihisto=1;
    else
        iunits=4;
        ihisto=0;
    end

    MODIS_varname2=remove_character(MODIS_varname,'_',' ');
    

%    nan_value=INFO.Vgroup(1).Vgroup(1).SDS(nvar).Attributes(2).Value;
    nan_value = fillval(var_nc);

    inan=find(modis_data_read==nan_value);
    modis_data_read(inan)=NaN;
        
%    add_offset=INFO.Vgroup(1).Vgroup(1).SDS(nvar).Attributes(iunits+2).Value;
    add_offset=var_nc.add_offset{1};
    modis_data_read=modis_data_read-add_offset; 
    %Not sure about the order to apply offsets and scale factors - offset
    %is zero for soil moisture anyway

%    scale_factor=INFO.Vgroup(1).Vgroup(1).SDS(nvar).Attributes(iunits+1).Value;
    scale_factor=var_nc.scale_factor{1};
%    if strcmp(scale_factor,'none')==1
%        scale_factor=1;
%    end
    modis_data_read=modis_data_read*scale_factor;

%    units_str=INFO.Vgroup(1).Vgroup(1).SDS(nvar).Attributes(iunits).Value;
    units_str=var_nc.units{1};
    
%% *** Make .timeseries3 for specific locations (keeping the NaNs) ***
    if exist('itimeseries_MODIS') & itimeseries_MODIS==1
        make_timeseries_SMOS
    end
    
    if ihisto==1    
        F = modis_data_read;
        F(isnan(F))=0;  %the NaN values within the histogram for each location
        %are kinda meaningless - think they should just mean they are zeros
        %(i.e. no pixels with those values). So set to zero here.
        totNpix = squeeze( sum(sum(F,1),2) ) ./ cf; %divide by cloud fraction to get the total sample pix (not just cloudy)
        igood=~isnan(totNpix); %when no. pixels is not NaN - i.e. when there were some pixels sampled
        inan=isnan(totNpix); %when no. pixels is not NaN
    else
        igood=~isnan(modis_data_read); %when no. pixels is not NaN
        inan=isnan(modis_data_read); %when no. pixels is not NaN
    end
    
    if exist('iaverage_modis') & iaverage_modis==1
        if imr~=1  %if averaging and not on the first file
                eval( ['Ndata = ' MODIS_varname '.Ndata + igood;'] ); %running total of no. days with some data
                
            if ihisto==1                
                eval( ['modis_data = ' MODIS_varname '.data + F;'] ); %add F here as we've removed all NaNs
                eval( ['modis_data_squared = ' MODIS_varname '.data_squared + F.^2;'] ); %add F here as we've removed all NaNs
                totNpix(inan)=0;
                eval( ['totNpix2 = ' MODIS_varname '.totNpix + totNpix;'] );
            else
                modis_data_read(inan)=0;  %inan only describes NaNs due to not being sampled that day 
                %(but histos have NaNs when have no pixels that were within that
                %range, which is why we use F).
                eval( ['modis_data = ' MODIS_varname '.data + modis_data_read;'] ); 
                eval( ['modis_data_squared = ' MODIS_varname '.data_squared + modis_data_read.^2;'] );                 
            end
        else
            Ndata = igood;
            modis_data_read(inan)=0;
            
            if ihisto==1                
                modis_data = F;
                modis_data_squared = F.^2;
                Ndata=igood;
                totNpix(inan)=0;
                totNpix2 = totNpix;
            else
                modis_data = modis_data_read; %this will be a running cumulative sum of the data
                modis_data_squared = modis_data_read.^2; %this will be a running total of data squared
            end
        end
    else
        modis_data = modis_data_read;
        modis_data_squared = 0;
        Ndata = 0;
        if ihisto==1
            totNpix2 = totNpix;
        end
    end

%% now put the data into the specific structure of name MODIS_varname
    eval([MODIS_varname '.data=modis_data;']);
    eval([MODIS_varname '.data_squared=modis_data_squared;']);    
    eval([MODIS_varname '.name2=MODIS_varname2;']);
    eval([MODIS_varname '.units_str=units_str;']);
    
    if exist('iaverage_modis') & iaverage_modis==1 & imr==nMOD_av
        eval([MODIS_varname '.modis_day_str=[num2str(nMOD_av) '' day TIME AVERAGE''];']);
        eval([MODIS_varname '.modis_year_str=modis_year_str;']);
    else
        eval([MODIS_varname '.modis_day_str=modis_day_str;']);
        eval([MODIS_varname '.modis_year_str=modis_year_str;']);
    end
    
    eval([MODIS_varname '.Ndata=Ndata;']);

    
    if ihisto==1        
        eval([MODIS_varname '.totNpix=totNpix2;']);
      
        if exist('iaverage_modis') & iaverage_modis==1 & imr==nMOD_av  %on the last loop of the time average divide by the number of datapoints used
            eval([MODIS_varname '.totNpix=' MODIS_varname '.totNpix ./ ' MODIS_varname '.Ndata;'])
            eval(['Ndata4D = repmat(' MODIS_varname '.Ndata,[1 1 size(modis_data,1) size(modis_data,2)]);']);
            Ndata4D=permute(Ndata4D,[3 4 1 2]);                                        
            
        end
    else
        eval(['Ndata4D=' MODIS_varname '.Ndata;']);
    end
    
    
    
    if exist('iaverage_modis') & iaverage_modis==1 & imr==nMOD_av  %on the last loop of the time average divide by the number of datapoints used
        eval(['dat_tot=' MODIS_varname '.data;']);
        eval(['dat_squared_tot=' MODIS_varname '.data_squared;']);   
        mean_dat = dat_tot ./ Ndata4D;
        std_dev = sqrt ( (dat_squared_tot -2*mean_dat.*dat_tot + Ndata4D.*mean_dat.^2) ./ (Ndata4D-1) );
        eval([MODIS_varname '.time_stdev = std_dev;']);
        eval([MODIS_varname '.data = mean_dat;']);                              
    end
    
   
end %for iread_modis


if exist('iaverage_modis')
        clear iaverage_modis
end

if exist('itimeseries_MODIS') & itimeseries_MODIS==1
       clear itimeseries_MODIS
end


disp('Done read MODIS');

