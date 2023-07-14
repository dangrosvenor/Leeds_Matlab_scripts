seaice = nc{'SEAICE'}(time,:,:);
landmask = nc{'LANDMASK'}(time,:,:);
terrain = nc{'HGT'}(time,:,:);

ice_shelf = landmask - seaice - terrain; %land is indicated whenever there is land or sea ice
                                                               %so if minus the seaice flag then whenever result is one have ice shelf
                                                               
save('c:/documents and settings/dan/my documents/WRF/ice_shelf_d03','ice_shelf');        

'done make ice shelf data'