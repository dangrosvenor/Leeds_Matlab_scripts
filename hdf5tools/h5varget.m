function data = h5varget ( varargin )
% H5VARGET:  reads data from an HDF5 file.
%
% DATA = H5VARGET(HFILE,HVAR) retrieves all of the data from the
% variable HVAR in the file HFILE.
%
% DATA = H5VARGET(HFILE,HVAR,START,COUNT) reads a contiguous subset
% of data from the variable HVAR in the file HFILE.  The contiguous 
% subset is a hyperslab defined by the index vectors START and
% COUNT, which are zero-based.
%
% DATA = H5VARGET(HFILE,HVAR,START,COUNT,STRIDE) reads a strided
% subset of data from the variable HVAR in the file HFILE.  The 
% strided subset will begin at the index START, have a length extent
% of COUNT along each dimension, but have an inter-element distance 
% given by STRIDE.  When not explicitly specified, STRIDE is 
% implicitly a vector of ones, i.e. contiguous data.
%
% DATA = H5VARGET(...,'enumerate',true) will transform the return
% values of an H5T_ENUM dataset as the char enumeration.  This 
% makes the return value a cell array of chars instead of the 
% numeric values.
%
% Example:
%   data = h5varget ( 'example.h5', '/g2/dset2.1' );  
%

%   Copyright 2007-2009 The MathWorks, Inc.

error(nargchk(2,7,nargin,'struct'));

args = parse_h5_args ( varargin{:} );

%
% Just use the defaults for now?
flags = 'H5F_ACC_RDONLY';
plist_id = 'H5P_DEFAULT';

file_id      = H5F.open ( args.filename, flags, plist_id );
dataset_id   = H5D.open(file_id,args.dataset);
datatype_id  = H5D.get_type(dataset_id);


dataspace_id = create_dataspace_id(args);
filespace_id = create_filespace_id(dataset_id,args);
memtype_id   = create_memtype_id(dataset_id,datatype_id);

data = H5D.read(dataset_id, memtype_id, dataspace_id, filespace_id, plist_id);

% Try to honor a fill value if the datatype is floating point.
if H5T.equal(datatype_id,'H5T_NATIVE_DOUBLE') || H5T.equal(datatype_id,'H5T_NATIVE_FLOAT')
	% Turn any fill value into NaN
	dcpl         = H5D.get_create_plist(dataset_id);
	if (H5P.fill_value_defined(dcpl) == H5ML.get_constant_value('H5D_FILL_VALUE_USER_DEFINED'))
		fillvalue    = H5P.get_fill_value(dcpl,datatype_id);
		data(data==fillvalue) = NaN;
	end
end


if ( H5T.get_class(datatype_id) == 8 ) 
	data = post_process_enums ( data, dataset_id, datatype_id, args );
end


%if H5T.equal ( datatype_id, H5T.copy('H5T_STD_REF_OBJ') )
%	data = post_process_references ( data );
%end

%
% Remove singleton dimensions.
data = squeeze(data);

H5T.close(datatype_id);
H5D.close(dataset_id);
H5F.close(file_id);


return










%-------------------------------------------------------------------------------
function memtype_id = create_memtype_id(dataset_id, datatype_id)

if ( H5T.get_class(datatype_id) == H5ML.get_constant_value('H5T_ENUM' ) )
	% If the datatype class is ENUM, then we need to recreate the memory datatype.
	memtype_id = H5T.copy(datatype_id);
	%memtype_id = 'H5ML_DEFAULT';
else
	% Otherwise, the default setting should do.
	memtype_id = 'H5ML_DEFAULT';
end

return



%===============================================================================
function data = post_process_enums ( data , dataset_id, datatype_id, args)

% Used to do this via preference, but that's only for backwards
% compatibility now.
if getpref('HDF5TOOLS', 'ENUMERATE', false) 
	args.enumerate = true;
end

if args.enumerate

	% Map from the values to a matlab index.
	idx = zeros(size(data));

	nmember = H5T.get_nmembers ( datatype_id );
	name = cell(nmember,1);
	for j = 0:nmember-1
		name{j+1} = H5T.get_member_name ( datatype_id, j );
		memb_value(j+1) = H5T.enum_valueof(datatype_id,name{j+1});
		idx(data==memb_value(j+1)) = j+1;
	end
	data = name(idx);

end

return

%===============================================================================
function dataspace_id = create_dataspace_id(args)

%
% Create the appropriate output dataspace.
if isempty(args.offset) && isempty(args.count) && isempty(args.stride)
	dataspace_id = 'H5S_ALL';
else
	dataspace_id = H5S.create_simple(length(args.offset),fliplr(args.count),fliplr(args.count));
end

return




%===============================================================================
%
% CREATE_FILESPACE_ID
%
% Create the dataspace that corresponds to the given selection.
function filespace = create_filespace_id(dataset_id,args)

%
% Create the appropriate mem dataspace
if isempty(args.offset) && isempty(args.count) && isempty(args.stride)
	filespace = 'H5S_ALL';
else
	%
	% define the memory hyperslab
	filespace = H5D.get_space(dataset_id);

	H5S.select_hyperslab(filespace, 'H5S_SELECT_SET', ...
	                     fliplr(args.offset), fliplr(args.stride), fliplr(args.count), ...
						 ones(1,length(args.offset)));

	% Is it valid?
	[ndims, dims] = H5S.get_simple_extent_dims(filespace);
	[bbsel_start,bbsel_finish] = H5S.get_select_bounds(filespace);
	if ( any((dims(:) - bbsel_finish(:)) < 0 ) ) 
		bb_selection_start_str = sprintf ( '%d ', bbsel_start );
		bb_selection_finish_str = sprintf ( '%d ', bbsel_finish );
		bb_selection_str = sprintf ( '[%s] to [%s]', deblank(bb_selection_start_str), ...
		                             deblank(bb_selection_finish_str) );
		bb_dataset_str = sprintf ( '%d ', dims );
		bb_dataset_str = sprintf ( '[%s]', deblank(bb_dataset_str) );
		error('HDF5TOOLS:H5VARGET:invalidFilespace', ...
		      'The bounding box on your selection (%s) exceeds the maximal extents of the dataset (%s).', ...
			  bb_selection_str, bb_dataset_str );
	end


end

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = post_process_references ( data )
% The result  of this is double precision.
%
% Leading dimension needs to be squeezed out?
power_of_16 = double(256 .^ (0:7));
data = power_of_16 * double(data);
return

