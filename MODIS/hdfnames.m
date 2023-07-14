[ ndata, nglob, stat]=hdfsd('fileinfo',sd_id);

for i=0:ndata-1
    sds_id = hdfsd('select',sd_id,i);
    
   %[dsname, dsndims, dsdims, dstype, dsatts, stat] =hdfsd('getinfo',sds_id)
              
	[dsname]= hdfsd('getinfo',sds_id);
    fprintf(1,'%d %s\n',i,dsname);
end