%read UM aerosol namelist and make an new one based on a constant scale
%factor, or if have new levels (or both)

scale_fac = 2.85; %superseeded if iheight_dep_scale_fac is set to one (sf set to 3.8 to go to 210 cm3 from 100 cm3)
iheight_dep_scale_fac=1; %E.g. once made a step change in the aerosol.to test the effect of a discontinuity.

%step-change
sf_ilev=[8 -1]; %Levels over which to apply sf_n. -1 indicates the end
sf_n=[1 0]; %Values to apply in the ranges specified by sf_ilev - should be of length length(sf_ilev)
    
%Linear decrease in height (of per kg value - so set flag to readin namelist as per kg not per m3)
x=[21:50];
sf_ilev=[20 x -1]; %Levels over which to apply sf_n. -1 indicates the end
val0=1.0;
ramp = val0 * (1 - (x - x(1)) / (x(end)-x(1)) );
sf_n=[1 ramp 1-val0]; %Values to apply in the ranges specified by sf_ilev - should be of length length(sf_ilev)
    

nlevs_new = 142;
nlevs_new = 70; %70 was the original no. levels

file_output = '/home/disk/eos1/d.grosvenor/UM/accum_100cm3_L142.nml';
file_output = '/home/disk/eos1/d.grosvenor/UM/accum_100cm3_L70_vertical_step_L8.nml';
file_output = '/home/disk/eos1/d.grosvenor/UM/accum_21cm3.nml';
file_output = '/home/disk/eos1/d.grosvenor/UM/accum_600cm3.nml';
file_output = '/home/disk/eos1/d.grosvenor/UM/accum_100cm3_Nd.nml';
file_output = '/home/disk/eos1/d.grosvenor/UM/accum_3.8e8_perkg_Annette_Aitken_Coarse_TAPER.nml';

%basefile = '/home/disk/eos1/d.grosvenor/UM/accum_21cm3.nml';
%basefile = '/home/disk/eos1/d.grosvenor/UM/accum_100cm3.nml';
%basefile = '/home/disk/eos1/d.grosvenor/UM/accum_210cm3.nml';
basefile = '/home/disk/eos1/d.grosvenor/UM/accum_3.8e8_perkg_Annette_Aitken_Coarse.nml';

basefile_zlevs = '/home/disk/eos1/d.grosvenor/UM/zlevs_orig_L70_40';
%The new file only gets used if nlevs ~= nlevs_new
file_zlevs_new = '/home/disk/eos1/d.grosvenor/UM/L142_40km_dz15m-200-1500m_zlevs';


fid=fopen(basefile,'rt');
dat= fscanf(fid,'%f',[1 inf]);
fclose(fid);

nlevs=70; %no. model levels that the base file is set up for
nmodes=10; %10 modes in total (mass and number).

clear dat2
for i=1:nmodes
    ii = (i-1)*(nlevs+1)+1; %index position for the start of the data for this mode - N.B. there is also an extra line at the top for each 
        %mode.
    dat2(:,i) = dat(:,ii:ii+nlevs);
end

if iheight_dep_scale_fac==0
    %Multiply all aerosol data by the scale factor
    %dat2(2,:) is the stash code, so leave this
    dat2(2:nlevs+1,:)=dat2(2:nlevs+1,:)*scale_fac;
else
    i1=1;
    for ilevs=1:length(sf_n)
        i0=i1+1; %start 1 ahead of the previous value
        i1=sf_ilev(ilevs)+1; %add one to the value specified since in sf_ilev since the first value is the stash code
        if i1==0 % started as -1, but will =0 now that have added one above
            i1=nlevs+1; %one extra because the first one is the stash code
        end
        dat2(i0:i1,:)=dat2(i0:i1,:)*sf_n(ilevs);
    end
end


if nlevs_new~=nlevs
   
    %read in original levels file - heights of levels
    fid=fopen(basefile_zlevs,'rt');
    tmp=fgetl(fid);
    dat_zlevs = textscan(fid,'%f %s %f %s %f %s %s\n');
    fclose(fid);
    
    N=dat_zlevs{1}(:);
    z=dat_zlevs{3}(:);
    
    %read in new levels file
    fid=fopen(file_zlevs_new,'rt');
    tmp=fgetl(fid);
    dat_zlevs = textscan(fid,'%f %s %f %s %f %s %s\n');
    fclose(fid);
    
    N2=dat_zlevs{1}(:);
    z2=dat_zlevs{3}(:);
    
    
    dat2_orig=dat2;
    clear dat2
    
    for i=1:nmodes
        dat2(1,i) = dat2_orig(1,i); %save the aerosol mode code=
        dat2(2:nlevs_new+1,i) = interp1(z,dat2_orig(2:nlevs+1,i),z2,'linear','extrap');                
    end
    
    
    
    
    
end


%% Now write out the result
fid=fopen(file_output,'wt');

for i=1:nmodes
    fprintf(fid,'%i\n',dat2(1,i));
    fprintf(fid,'%.5e\n',dat2(2:end,i));    

end
fclose(fid);


