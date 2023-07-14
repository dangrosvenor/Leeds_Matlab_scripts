data_type='wind';
%data_type='potemp';
data_type='vert wind';

XY_pot_cross_data.X_cross = timesTH(1).t(1:pend);
XY_pot_cross_data.Y_cross = zz(1).z; 

switch data_type
    case 'wind'
        U_cross_6th_Jan_6UTC = pdat(1).p(:,1:pend);
    case 'potemp'
        pot_cross_6th_Jan_6UTC = pdat(1).p(:,1:pend);
    case 'vert wind'
        vertwind_cross_6th_Jan_6UTC = pdat(1).p(:,1:pend);

end

savedir = 'Y:\WRF\ecmwf_ml_0.5_nudging\figures_for_paper_Aug2010\';
savename = [savedir 'pot_slice_streamline_LAT=67.5_z0=950_str_succ=129_0-5km.mat'];
isave=1;
if isave==1
    save(savename,'XY_pot_cross_data','U_cross_6th_Jan_6UTC','pot_cross_6th_Jan_6UTC',...
        'lon_slice','d_line');
end