funciton new_data=map_data_onto_cticks(data,xvals,cvals)
%funciton new_data=map_data_onto_cticks(data,xvals,cvals)
%maps data onto new colorscale for supplied tick mark values
% xvals vary from 0 to 1 and specify the distance along the colorbar
% cvals are the corresponding data values
% linearly interpolates between the values

new_data=interp1(cvals,xvals,data(:));