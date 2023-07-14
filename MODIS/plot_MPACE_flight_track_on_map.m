%plot the MODIS aircraft flight track

figure
%m_proj('lambert','lon',[-180 -135],'lat',[68 74]);
%m_proj('lambert','lon',[-158 -147],'lat',[70 71.5]);
m_proj('lambert','lon',[-151 -149],'lat',[70.3 70.6]);

m_plot(mpace_dat(col_lon,:),mpace_dat(col_lat,:),'b-');
m_grid('box','fancy','tickdir','in');
hcoast=m_coast('line','linewidth',2,'color','k');

%savedir='C:\Users\Dan\Documents\plots\';
savedir='C:\Users\Dan\Documents\logbook\Antarctica\Flights and instruments\plots\';
savename=[savedir 'Flight_track ' flight_no];