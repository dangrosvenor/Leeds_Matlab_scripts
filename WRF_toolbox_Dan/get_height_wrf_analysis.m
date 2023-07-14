function height = get_height_wrf_analysis(nc,ilat,ilon,iloc)

                if prod(size(nc{'GHT'})) > 0
                    height = nc{'GHT'}(1,:,ilat(iloc),ilon(iloc));
                    incep=1;
                else
                    incep=0;

                    pres=nc{'PRES'}(1,:,ilat(iloc),ilon(iloc)); %BUT should the first pressure and temperature be counted - know
                    % that the pressure is the surface pressure - what does that mean for the temperature - skintemp?
                    pres=pres(2:end);   %don't want the first level - at least for ecmwf ml data
                    psfc=nc{'PSFC'}(:); psfc=psfc(ilat(iloc),ilon(iloc));
                    pmsl=nc{'PMSL'}(:); pmsl=pmsl(ilat(iloc),ilon(iloc));
                    temp=nc{'TT'}(1,:,ilat(iloc),ilon(iloc));
                    temp=temp(2:end);

                    %Height of the terrain as input from the analysis (underlying) model (coarse resolution interpolated onto wrf grid).
                    % (PSFC - this is the pressure at the height of the analysis terrain
                    % and not the pressure at the height of the high resolution terrain
                    %as taken from the geogrid program and given in HGT_M - tested this by checking the PSFC field in met_em
                    %files and wrfout and they are consistent with this).
                    %*** Actually, looking at the first level of PRES and comparing it to PSFC it seems the two are very close with both positive
                    %and negative differences. the PSFC lines look less smooth. PRES seems to follow the terrain more closely in terms of this
                    %smoothness anf the peaks and troughs correspond more closely - so this may be the pressure to choose to
                    %correspond to the height of SOILHGT
                    %NOTE that SOILHGT can be negative
                    %difference between choosing PSFC and PRES level one is about 80-100 m in the troposphere

                    %Here the hydrostatic equation is solved starting at the surface pressure, which is known to be at the height
                    %of the analysis terrain. The temperature vs pressure profile is described by PRES and TT and the required temperature
                    %at a given pressure (e.g. PSFC to start with) is found by interpolation in the hydrostatic2 function

                    %NOTE, might want to use this for finding pressure corrections for NCEP analysis files
                    %                     soilhgt = nc{'GHT'}(1,:,ilat(iloc),ilon(iloc)); %height of the second (=1000 hPa) pressure level
                    %                     soilhgt = soilhgt(2);
                    %                     [Y,I]=sort(pres,'descend'); %sort the data in order to take into account pres(2)(=surface) being lower than first level of 1000 hPa
                    %                     pres=pres(I);
                    %                     temp=temp(I);


                    if prod(size(nc{'SOILHGT'})) > 0
                        soilhgt=nc{'SOILHGT'}(:); soilhgt=soilhgt(ilat(iloc),ilon(iloc));  %height of PSFC level
                        pres2=pres;
                        temp2=temp;

                    else
                        pres2 = [pmsl pres];
                        skinT = temp(1)*(pmsl/pres(1))^0.286; %temperature assuming a constant potential temperature between msl and the surface point
                        %skinT = nc{'SKINTEMP'}(:); %using skin temperature as sea level temp - could be inaccurate?
                        %skinT = skinT(ilat(iloc),ilon(iloc));
                        %    skinT=skinT+5; %to test temperature sensitivity - does not seem that senstitive - 5 degree increase led to only 0.14 hPa increase in
                        %in estimated pressure of a height below soilhgt (tried 91 m wheras soilhgt was 111m - difference was even less for lower altitudes).
                        temp2 = [skinT temp];
                        %                soilhgt=nc{'SOILHGT'}(:); soilhgt=soilhgt(ilat(iloc),ilon(iloc));  %height of PSFC level
                        soilhgt=0;
                    end




                    PSPAN=[pres2(1) pres2(end)]; %pressure range for integration
                    [P,hp] = ODE45(@hydrostatic2,PSPAN,soilhgt,[],pres2,temp2); %enter the initial height for the given the first value of PSPAN
                    %note that some of the pressure levels in the PRES array will be below the

                    height = interp1(P,hp,pres);

                    %                     pres=nc{'PRES'}(time,:,ilat(iloc),ilon(iloc))/100;
                    %                     height = find_height_from_p_ant_d03(pres);

                end