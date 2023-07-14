function SummarizeSatelliteGeometry(satGEOM_struct,SatID,origin_llh,...
                                      tsince_min,plotTitle)

dtr=pi/180;
min_per_day=60*24;

sat_elev =satGEOM_struct.sat_elev;
sat_llh  =satGEOM_struct.sat_llh;
sat_rnge =satGEOM_struct.sat_rnge;
sat_rdot =satGEOM_struct.sat_rdot;
sat_elev =satGEOM_struct.sat_elev;
sat_phi  =satGEOM_struct.sat_phi;
satp_tcs =satGEOM_struct.satp_tcs;
satp_llh =satGEOM_struct.satp_llh;
xyzp     =satGEOM_struct.xyzp;
%uk_xyzp=satGEOM_struct.
rngp     =satGEOM_struct.rngp;
thetap   =satGEOM_struct.thetap;
phip     =satGEOM_struct.phip;
vp       =satGEOM_struct.vp;
vk       =satGEOM_struct.vk;
%s=satGEOM_struct.
thetaB  =satGEOM_struct.thetaB;
phiB    =satGEOM_struct.phiB;
cosBP=satGEOM_struct.cosBP;

sat_utsec=(tsince_min-tsince_min(1))*60;
UThrs=tsince_min/60;

fprintf('Summary data from file %s \n',plotTitle)
fprintf('Station %s: Lon=%6.4f Lat %6.4f Alt=%6.2f m \n\n',...
    'Kwaj',origin_llh(2)/dtr,origin_llh(1)/dtr,origin_llh(3))
fprintf(' YEAR MO  DAY UThrs \n')
fprintf('%5i %2i %4i %5.2f \n',SatID.year,SatID.mon,SatID.day,UThrs(1));
fprintf('Pass Duration %6.2f min %6.2f degrees \n\n',...
    (UThrs(end)-UThrs(1))/60,max(sat_elev/dtr))

Display=input('Input 1 for station geometry summary, else CR \n');
if Display
    figure
    plot(sat_llh(2,:)/dtr,sat_llh(1,:)/dtr,'b')
    hold on
    plot(satp_llh(2,:)/dtr,satp_llh(1,:)/dtr,'r')
    hold on
    plot(origin_llh(2)/dtr,origin_llh(1)/dtr,'mp')
    grid on
    legend('Satellite','300 km intercept')
    xlabel('Longitude--deg')
    ylabel('Latitude--deg')
    title(plotTitle)
    bold_fig
    
    figure
    subplot(2,1,1)
    plot(sat_utsec/3600,sat_rnge/1000,'r')
    grid on
    title(plotTitle)
    ylabel('Range--km')
    subplot(2,1,2)
    plot(sat_utsec/3600,sat_rdot/1000,'b')
    grid on
    title(plotTitle)
    xlabel('UT-hrs')
    ylabel('Rangerate--km/s')
    bold_fig
    
    figure
    plot(sat_utsec/3600,sat_elev/dtr,'r')
    hold on
    plot(sat_utsec/3600,sat_phi/dtr,'b')
    grid on
    legend('Elevation-deg','Azimuth--deg')
    title(plotTitle)
    xlabel('UT-hrs')
    ylabel('angle--deg')
    bold_fig
end


Display=input('Input 1 for penetration point geometry summary, else CR \n');
if Display
    figure
    subplot(2,1,1)
    plot(satp_llh(2,:)/dtr,satp_llh(1,:)/dtr,'b')
    hold on
    plot(satp_llh(2,1)/dtr,satp_llh(1,1)/dtr,'b.')
    hold on
    plot(origin_llh(2)/dtr,origin_llh(1)/dtr,'mp')
    ylabel('Latitude--deg')
    xlabel('Longitude--deg')
    grid on
    title(plotTitle)
    subplot(2,1,2)
    plot(sat_utsec/3600,thetap/dtr,'r')
    hold on
    plot(sat_utsec/3600,phip/dtr,'b')
    grid on
    title(plotTitle)
    xlabel('UT--hrs')
    ylabel('degrees')
    legend('\theta','\phi')
    bold_fig
end


Display=input('Input 1 to display penetration point velocity, else CR \n');
if Display
    figure
    subplot(3,1,1)
    plot(sat_utsec/3600,vp(1,:)/1000,'r')
    grid on
    ylabel('vpx--kms')
    title(plotTitle)
    subplot(3,1,2)
    plot(sat_utsec/3600,vp(2,:)/1000,'r')
    grid on
    ylabel('vy--kms')
    subplot(3,1,3)
    plot(sat_utsec/3600,vp(3,:)/1000,'r')
    grid on
    ylabel('vpz--kms')
    xlabel('UT--hrs')
    bold_fig
end

Display=input('Input 1 to display vk, else CR \n');
if Display
    figure
    subplot(3,1,1)
    plot(satp_llh(2,:)/dtr,satp_llh(1,:)/dtr,'b')
    hold on
    plot(satp_llh(2,1)/dtr,satp_llh(1,1)/dtr,'b.')
    hold on
    plot(origin_llh(2)/dtr,origin_llh(1)/dtr,'mp')
    grid on
    legend('Penetration Point')
    ylabel('Latitude--deg')
    title(plotTitle)
    subplot(3,1,2)
    plot(sat_utsec/3600,vk(1,:)/1000,'r')
    grid on
    ylabel('vk-east--kms')
    subplot(3,1,3)
    plot(sat_utsec/3600,vk(2,:)/1000,'r')
    grid on
    ylabel('vk-south--kms')
    xlabel('UT--hrs')
    bold_fig
end

Display=input('Input 1 to display B field geometry, else CR \n');
if Display
    
    if  exist('Bz100','file')
    load('Bz300')  %Generated by DemoIGRF
    figure
    imagesc(zlon_deg,zlat_deg,dip'/dtr)
    axis xy
    hold on
    plot(satp_llh(2,:)/dtr,satp_llh(1,:)/dtr,'mp')
    hold on
    plot(origin_llh(2)/dtr, origin_llh(1)/dtr,'wp')
    colorbar
    title('Polar angle from vertical reference plane')
    xlabel('Longitude--deg')
    ylabel('Latitude--deg')
    LINESPEC='w';
    hold on
    [tempdat,h]=contour(zlon_deg,zlat_deg,dip'/dtr,20,LINESPEC);
    end
    
    figure
    plot(satp_llh(2,:)/dtr,satp_llh(1,:)/dtr,'b')
    hold on
    plot(satp_llh(2,1)/dtr,satp_llh(1,1)/dtr,'b.')
    hold on
    plot(origin_llh(2)/dtr,origin_llh(1)/dtr,'mp')
    ylabel('Latitude--deg')
    grid on
    title(plotTitle)
    
    figure
    plot(sat_utsec/3600,thetaB/dtr,'r')
    hold on
    plot(sat_utsec/3600,phiB/dtr,'b')
    grid on
    axis([floor(10*sat_utsec(1))/36000 ceil(10*sat_utsec(end))/36000 -200 100])
    legend('\theta-B','\phi-B')
    title(plotTitle)
    
    figure
    plot(sat_utsec/3600,cosBP,'r')
    grid on
    axis([floor(10*sat_utsec(1))/36000 ceil(10*sat_utsec(end))/36000 -1 1])
    ylabel('cosBP')
    xlabel('UT--hrs')
    title(plotTitle)
end
return