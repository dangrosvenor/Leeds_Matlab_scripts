function data=read_h4_var(SD_id,nvar)

sds_id = hdfsd('select',SD_id, nvar); %get the variable ID for variable nvar
[name,rank,dimsizes,data_type,nattrs,status] = hdfsd('getinfo',sds_id); %get the info for this variable
start=zeros(size(dimsizes));
stride=ones(size(dimsizes));
edge=dimsizes;
[data,status] = hdfsd('readdata',sds_id,start,stride,edge); %get the data
data=double(data);
