%read the cat2 data with ceilometer data included

ceil_file = ['/home/disk/eos1/d.grosvenor/modis_work/Irshad_data/ML_SL_cloud_cat2.mat'];
load(ceil_file);
% 
% Last column 40-49 is Ceilometer data. And columns before that are your
% cat 2 data
% 
% Column 46 = 1st cloud layer
% 
% Column 47 = 2nd cloud layer
% 
% Column 48 = 3rd cloud layer
% 
% Column 49 = 4th cloud layer

%values are set to zero by default - i.e. when there is not cloud base
%height for a given layer.

%now match to the available MODIS data - note that this ceilometer data is
%only available in this file for the cat 2 CTT>265K matches and not all of the MODIS
%data. Irshad should be able to generate matches for all data if required.

%note this data is arrange year, month, day, which is what datenum requries
%But note that the data I sent Irshad the date is arranged year, day, month

Matlab_time_ceil = datenum(Res_Dan_800(:,1),Res_Dan_800(:,2),Res_Dan_800(:,3),Res_Dan_800(:,4),Res_Dan_800(:,5),0);

ceil_layer01 = NaN*ones(size(Date_Time_Swath.timeseries3));
ceil_layer02 = NaN*ones(size(Date_Time_Swath.timeseries3));
ceil_layer03 = NaN*ones(size(Date_Time_Swath.timeseries3));
ceil_layer04 = NaN*ones(size(Date_Time_Swath.timeseries3));
for i=1:length(Matlab_time_ceil)
    itime = find(abs(Matlab_time_ceil(i)-Date_Time_Swath.timeseries3)<0.1/60/24); %match to 0.1 mins in case of rounding errors
    if length(itime>0)
        ceil_layer01(itime) = Res_Dan_800(i,46);
        ceil_layer02(itime) = Res_Dan_800(i,47);
        ceil_layer03(itime) = Res_Dan_800(i,48);
        ceil_layer04(itime) = Res_Dan_800(i,49);        
    end        
end

ceil_all_layers = cat(4,ceil_layer01,ceil_layer02,ceil_layer03,ceil_layer04);
ceil_layer_maxH = max( ceil_all_layers , [],4);
ceil_nlayers_all = zeros(size(ceil_all_layers));
ceil_nlayers_all(ceil_all_layers>1e-5) = 1; %set to one when have a layer reported
ceil_nlayers = sum(ceil_nlayers_all,4); %sum the zeros and ones to get the number of layers




