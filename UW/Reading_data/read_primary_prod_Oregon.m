%Files only contain one field - npp (Net Primary Productivity)
filename_h5='/home/disk/eos8/d.grosvenor/primary_prod_Oregon/vgpm.2014335.hdf'; %This is a monthly mean. Presumably 335 rfers to day 335
 %, which would be 1st Dec, 2014. So this is for December.
INFO = hdfinfo(filename_h5);
SD_id = hdfsd('start',filename_h5,'read'); %open the file

%[npp,dimsizes_npp,sampling_along,sampling_across]=get_hdf_data_dan('npp',SD_id,INFO);
var_name = 'npp';

nvar = hdfsd('nametoindex',SD_id,var_name);
sds_id = hdfsd('select',SD_id, nvar);
[name,rank,dimsizes,data_type,nattrs,status] = hdfsd('getinfo',sds_id);
start=zeros([1 length(dimsizes)]);
stride=ones([1 length(dimsizes)]);
edge=dimsizes;
[npp,status] = hdfsd('readdata',sds_id,start,stride,edge);

npp = double(npp);

qpcolor(flipdim(npp',1)); shading flat; caxis([0 300]); colorbar;