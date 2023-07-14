function [ds_data,status]=gethdfdat(sd_id,idat)

sds_id = hdfsd('select',sd_id,idat);

[dsname, ds_ndims, ds_dims, dstype, dsatts, stat] ...
    = hdfsd('getinfo',sds_id);

ds_start = zeros(1,ds_ndims);
ds_stride=[];
ds_edges=ds_dims;

[ds_data, status] =hdfsd('readdata',sds_id,ds_start,ds_stride,ds_edges);