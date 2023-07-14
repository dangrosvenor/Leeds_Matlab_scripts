case_exec = 'stream';


switch case_exec
    case 'grid'
        %load on wrfout file for d03 - done in loadWRFvar now
        x_grid = ([1:size(lat2d(1).var,2)]-1) * distlatlon(lat2d.var(1,1),lon2d.var(1,1),lat2d.var(1,2),lon2d.var(1,2));
        y_grid = ([1:size(lat2d(1).var,1)]-1) * distlatlon(lat2d.var(1,1),lon2d.var(1,1),lat2d.var(2,1),lon2d.var(2,1));

    case 'stream'
        %then load the regular grid file
        
        zmax=4500; %max height of data to load in up to (takes a long time to load)
        
        iload=1; %in case have already loaded them set iload to zero
        if iload==1
            
            imanual_select_load_case=1;
            dire(1).dir = 'Y:\WRF\';
            rundir(1).dir='ecmwf_ml_0.5_nudging';
            fileWRF(1).file=['regular_zgrid_output_time11.nc'];   %             
            load_WRF_vars

            zfine=nc{'ZFINE'}(:);            
            izmax=findheight(zfine,zmax);
            
            ufine=nc{'UFINE'}(1,1:izmax,:,:); %only loading in the lower part of the domain as takes long time
            vfine=nc{'VFINE'}(1,1:izmax,:,:);
            wfine=nc{'WFINE'}(1,1:izmax,:,:);
            
            ufine=permute(ufine,[2 3 1]); %need to check that have the right order for x and y
            vfine=permute(vfine,[2 3 1]);
            wfine=permute(wfine,[2 3 1]);

        end
        
        %load the WRF file as load_WRF_vars calculates x_grid
        imanual_select_load_case=1;
        dire(1).dir = 'Y:\WRF\';
        rundir(1).dir='ecmwf_ml_0.5_nudging'; 
        fileWRF(1).file=['d03'];   %usual ECMWF run (24th April, 2009)     
        load_WRF_vars

        
        [xstream,ystream,zstream]=meshgrid(x_grid,y_grid,zfine(1:izmax));

        xstream=xstream*1000;
        ystream=ystream*1000; %convert to m

        %NOTE - the streamlines are calculated and drawn in streamlines_threeD_draw.m

%        S=stream3(xstream,ystream,zstream,ufine,vfine,wfine,200e3,400e3,1.25e3);
        %here the last 3 are the starting x,y,z positions for the streamline
        %looks good from this position - need to decide where is best to start the streamline

%        plot(S{1}(:,1)/1e3,S{1}(:,2)/1e3,'g'); %to plot the x,y co-ords

end

threeD_streamline_terrain_calc %works out the actual terrain height available from
%the fine grid - i.e. the lowest level at which data is present (not NaN) in ufine, etc.
%and calls it terr_fine

inan=find(ufine<-1e9); ufine2=ufine; ufine2(inan)=NaN; %stream3 just stops the streamlines when it reaches a NaN
inan=find(ufine<-1e9); vfine2=vfine; vfine2(inan)=NaN;
inan=find(ufine<-1e9); wfine2=wfine; wfine2(inan)=NaN;

dz=100;
diz=round(dz/(zfine(2)-zfine(1)));
inds=[diz:diz:size(ufine,3)];
ufine3=ufine2(:,:,inds);
vfine3=vfine2(:,:,inds);
wfine3=wfine2(:,:,inds);
zstream3=zstream(:,:,inds);
ystream3=ystream(:,:,inds);
xstream3=xstream(:,:,inds);

disp('Finished 3D streamline calc');
