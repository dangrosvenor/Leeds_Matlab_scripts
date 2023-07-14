flags = 'H5F_ACC_RDWR'; 
plist_id = 'H5P_DEFAULT'; 
varname = '/DSC1'; 
file_id     = H5F.open('foo.h5',flags,plist_id); 
datasetId  = H5D.open(file_id,varname); 
datatype_id = H5D.get_type(datasetId); 
dcpl = H5D.get_create_plist(datasetId); 
fileSpaceId = H5D.get_space(datasetId); 
% select the hyperslab to start at element #5 and extend for 
% 10 elements 
H5S.select_hyperslab(fileSpaceId, 'H5S_SELECT_SET',5,1,10,1); 
 
 % Now extend the dataset. 
 H5D.set_extent(datasetId,15); 
% H5D.write(datasetId, 'H5ML_DEFAULT', 'H5S_ALL', fileSpaceId, plist_id, [6:15]);
H5D.close(datasetId);
H5F.close(file_id);
