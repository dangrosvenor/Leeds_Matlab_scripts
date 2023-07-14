function netcdf_add_variable(file_netcdf,varstruct,nc)

varstruct2.Name = varstruct.Name;
varstruct2.Nctype = varstruct.Nctype;
varstruct2.Dimension = varstruct.Dimension;
nc_addvar(file_netcdf,varstruct2);

nc_attput(file_netcdf,varstruct.Name,'FieldType',varstruct.FieldType);
nc_attput(file_netcdf,varstruct.Name,'MemoryOrder',varstruct.MemoryOrder);
nc_attput(file_netcdf,varstruct.Name,'units',varstruct.units);
nc_attput(file_netcdf,varstruct.Name,'description',varstruct.description);
nc_attput(file_netcdf,varstruct.Name,'stagger',varstruct.stagger);


data = nc{varstruct.Name}(:);
nc_varput(file_netcdf,varstruct.Name,data);