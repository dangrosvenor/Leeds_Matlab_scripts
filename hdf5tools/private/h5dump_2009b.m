function h5dump_2009b(filename,objName)
% H5DUMP:  prints hdf5 metadata
%
%     H5DUMP(FILENAME) prints the entire file's metadata.
%
%     H5DUMP(FILENAME,OBJ) prints the metadata for the named object
%     object.  If obj is a group, all objects below the group will
%     be described.

% Copyright 2007-2009 The MathWorks, Inc.

global indent_char;

indent_char = '   ';
close_it = false;

idx_type = 'H5_INDEX_NAME';
order = 'H5_ITER_INC';

file_id = H5F.open(filename, 'H5F_ACC_RDONLY', 'H5P_DEFAULT');
if nargin == 1
	gid = H5G.open(file_id,'/');
	dump_group(gid);
	H5G.close(gid);
	H5F.close(file_id);
else
	oid = H5O.open(file_id,objName,'H5P_DEFAULT');
	info = H5O.get_info(oid);
	switch ( info.type )
		case H5ML.get_constant_value('H5O_TYPE_GROUP')
			dump_attributes(oid,1);
			dump_group(oid);
		case H5ML.get_constant_value('H5O_TYPE_DATASET')
			dump_dataset(oid);
		otherwise
			fprintf('unknown object type %d', info.type );
	end
	H5O.close(oid);
end




%-----------------------------------------------------------------------------
function [status opdata_out]= H5Literator_func(gid,name,opdata)
opdata_out = opdata;
status = 0;
try
	obj_id = H5O.open(gid,name,'H5P_DEFAULT');
catch me
	if strcmp(me.identifier,'MATLAB:hdf5lib2:H5Oopen:failure')
		return
	else
		rethrow(me);
	end
end
info = H5O.get_info(obj_id);

switch ( info.type )
	case H5ML.get_constant_value('H5O_TYPE_GROUP')
		dump_group(obj_id);
	case H5ML.get_constant_value('H5O_TYPE_DATASET')
		dump_dataset(obj_id);
	otherwise
		fprintf('unknown object type %d', info.type );
end
H5O.close(obj_id);

status = 0;
return


%-----------------------------------------------------------------------------
function dump_group(gid)

gname = H5I.get_name(gid);
fprintf ('\nGroup ''%s'' {\n', gname);
dump_attributes(gid,1);
H5L.iterate(gid, 'H5_INDEX_NAME', 'H5_ITER_INC', 0, @H5Literator_func, []);
fprintf ('} End Group ''%s''\n', gname);

%-----------------------------------------------------------------------------
function dump_dataset(dset_id)

global indent_char;

dset_name = H5I.get_name(dset_id);
fprintf ('%sDataset ''%s'' {\n', indent_char, dset_name);

dtype = H5D.get_type(dset_id);
dump_datatype(dtype,2);
H5T.close(dtype);

space_id = H5D.get_space(dset_id);
dump_dataspace(space_id,2);
H5S.close(space_id);

dump_layout(dset_id,2);
dump_attributes(dset_id,2);

fprintf ('%s} End Dataset ''%s''\n\n', indent_char, dset_name);

%--------------------------------------------------------------------------
function dump_attributes(loc_id,indent_level)
total_attrs = H5A.get_num_attrs(loc_id);
for j = 0:total_attrs-1
    attr_id = H5A.open_idx(loc_id,j);
    dump_attribute(attr_id,indent_level);
    H5A.close(attr_id);
end
return

%--------------------------------------------------------------------------
function dump_attribute(attr_id,indent_level)
global indent_char;
indent = repmat(indent_char, 1,indent_level);
fprintf('%sAttribute: ''%s'' {\n', indent, H5A.get_name(attr_id) );

dtype = H5A.get_type(attr_id);
dump_datatype(dtype,indent_level+1);
H5T.close(dtype);

space_id = H5A.get_space(attr_id);
dump_dataspace(space_id,indent_level+1);
H5S.close(space_id);

fprintf('%s}\n',indent);
return

%--------------------------------------------------------------------------
function dump_datatype(dtype, indent_level)

global indent_char;
indent = repmat(indent_char, 1,indent_level);

%dclass = H5T.get_class(dtype);
%sz = H5T.get_size(dtype);

dtypeStr = get_datatype_string(dtype,indent_level);

fprintf ('%sDatatype:  %s\n', indent, dtypeStr);

return

%--------------------------------------------------------------------------
function dump_dataspace(space_id,indent_level)

global indent_char;
indent = repmat(indent_char, 1,indent_level);

[ndims,dims,maxdims] = H5S.get_simple_extent_dims(space_id);
if ndims == 0
    fprintf ('%sDataspace scalar\n', indent);
else
    str = sprintf('%.0f ', flipud(dims(:)));
    str = deblank(str);

    %unlimited = H5ML.get_constant_value('H5S_UNLIMITED');
    strmax = '';
    for j = ndims:-1:1
        if maxdims(j) == 18446744073709551616 % what the hey??
            strmax = [strmax 'H5S_UNLIMITED ']; %#ok<AGROW>
        else
            strmax = [strmax num2str(maxdims(j)) ' ']; %#ok<AGROW>
        end
    end
    strmax = deblank(strmax);

    fprintf ('%sDataspace simple {[%s] / [%s]}\n', indent, str, strmax);
end

return

%--------------------------------------------------------------------------
function dtypeStr = get_datatype_string(dtype,indent_level)



% Check vlen strings first.
if H5T.is_variable_str(dtype)
    dtypeStr = get_vlen_string(dtype,indent_level);
	return
end

class_id = H5T.get_class(dtype);

switch ( class_id )
    case H5ML.get_constant_value('H5T_BITFIELD')
        if H5T.equal(dtype,'H5T_STD_B8BE')
            dtypeStr = 'H5T_STD_B8BE';
        elseif H5T.equal(dtype,'H5T_STD_B8LE')
            dtypeStr = 'H5T_STD_B8LE';
        elseif H5T.equal(dtype,'H5T_STD_B16LE')
            dtypeStr = 'H5T_STD_B16LE';
        elseif H5T.equal(dtype,'H5T_STD_B16BE')
            dtypeStr = 'H5T_STD_B16BE';
        elseif H5T.equal(dtype,'H5T_STD_B32LE')
            dtypeStr = 'H5T_STD_B32LE';
        elseif H5T.equal(dtype,'H5T_STD_B32BE')
            dtypeStr = 'H5T_STD_B32BE';
        elseif H5T.equal(dtype,'H5T_STD_B64LE')
            dtypeStr = 'H5T_STD_B64LE';
        elseif H5T.equal(dtype,'H5T_STD_B64BE')
            dtypeStr = 'H5T_STD_B64BE';
        else
            error('HDF5TOOLS:h5dump:unrecognizedBitFieldType', ...
                  'Encountered an unrecognized bitfield type.' );
        end
    case H5ML.get_constant_value('H5T_INTEGER')
        if H5T.equal(dtype,'H5T_STD_I8BE')
            dtypeStr = 'H5T_STD_I8BE';
        elseif H5T.equal(dtype,'H5T_STD_I8LE')
            dtypeStr = 'H5T_STD_I8LE';
        elseif H5T.equal(dtype,'H5T_STD_U8BE')
            dtypeStr = 'H5T_STD_U8BE';
        elseif H5T.equal(dtype,'H5T_STD_U8LE')
            dtypeStr = 'H5T_STD_U8LE';
        elseif H5T.equal(dtype,'H5T_STD_I16BE')
            dtypeStr = 'H5T_STD_I16BE';
        elseif H5T.equal(dtype,'H5T_STD_I16LE')
            dtypeStr = 'H5T_STD_I16LE';
        elseif H5T.equal(dtype,'H5T_STD_U16BE')
            dtypeStr = 'H5T_STD_U16BE';
        elseif H5T.equal(dtype,'H5T_STD_U16LE')
            dtypeStr = 'H5T_STD_U16LE';
        elseif H5T.equal(dtype,'H5T_STD_I32BE')
            dtypeStr = 'H5T_STD_I32BE';
        elseif H5T.equal(dtype,'H5T_STD_I32LE')
            dtypeStr = 'H5T_STD_I32LE';
        elseif H5T.equal(dtype,'H5T_STD_U32BE')
            dtypeStr = 'H5T_STD_U32BE';
        elseif H5T.equal(dtype,'H5T_STD_U32LE')
            dtypeStr = 'H5T_STD_U32LE';
        elseif H5T.equal(dtype,'H5T_STD_I64BE')
            dtypeStr = 'H5T_STD_I64BE';
        elseif H5T.equal(dtype,'H5T_STD_I64LE')
            dtypeStr = 'H5T_STD_I64LE';
        elseif H5T.equal(dtype,'H5T_STD_U64BE')
            dtypeStr = 'H5T_STD_U64BE';
        elseif H5T.equal(dtype,'H5T_STD_U64LE')
            dtypeStr = 'H5T_STD_U64LE';
        else
            error('HDF5TOOLS:h5dump:unrecognizedIntegerType', ...
                  'Encountered an unrecognized integer type.' );
        end

    case H5ML.get_constant_value('H5T_FLOAT')
        if H5T.equal(dtype,'H5T_IEEE_F32BE')
            dtypeStr = 'H5T_IEEE_F32BE';
        elseif H5T.equal(dtype,'H5T_IEEE_F32LE')
            dtypeStr = 'H5T_IEEE_F32LE';
        elseif H5T.equal(dtype,'H5T_IEEE_F64BE')
            dtypeStr = 'H5T_IEEE_F64BE';
        elseif H5T.equal(dtype,'H5T_IEEE_F64LE')
            dtypeStr = 'H5T_IEEE_F64LE';
        else
            error('HDF5TOOLS:h5dump:unrecognizedFloatType', ...
                  'Encountered an unrecognized floating point type.' );
        end

    case H5ML.get_constant_value('H5T_OPAQUE')
        dtypeStr = get_opaque_string(dtype,indent_level);

    case H5ML.get_constant_value('H5T_REFERENCE')
        dtypeStr = 'H5T_REFERENCE';

    case H5ML.get_constant_value('H5T_COMPOUND')
        dtypeStr = get_compound_string(dtype,indent_level);


    case H5ML.get_constant_value('H5T_ENUM')
        dtypeStr = get_enum_string(dtype,indent_level);


    case H5ML.get_constant_value('H5T_ARRAY')
        dtypeStr = get_array_string(dtype,indent_level);

    case H5ML.get_constant_value('H5T_VLEN')
        dtypeStr = dump_vlen_datatype(dtype,indent_level+1);

    case H5ML.get_constant_value('H5T_STRING')
        dtypeStr = 'H5T_STRING';

    otherwise
        error('HDF5TOOLS:h5dump:unrecognizedClass', ...
              'Encountered an unrecognized datatype class %d.', ...
              class_id);
end


return

%--------------------------------------------------------------------------
function dump_layout(dset_id,indent_level)

global indent_char;
indent = repmat(indent_char, 1,indent_level);

dcpl = H5D.get_create_plist(dset_id);

layout = H5P.get_layout(dcpl);
switch ( layout )
    case H5ML.get_constant_value('H5D_COMPACT')
		layout_string = sprintf('Layout:  Compact');

    case H5ML.get_constant_value('H5D_CONTIGUOUS')
		layout_string = sprintf('Layout:  Contiguous');

    case H5ML.get_constant_value('H5D_CHUNKED')
		layout_string = sprintf('Layout:  Chunked');
		[rank,dims] = H5P.get_chunk(dcpl);
		str = '';
    	for j = rank:-1:1
	        str = [str num2str(dims(j)) ' ']; %#ok<AGROW>
	    end
		str = sprintf('{[%s]}', deblank(str));
		layout_string = sprintf ( '%s %s', layout_string, str);

end

fprintf ('%s%s\n', indent, layout_string);
H5P.close(dcpl);

return


