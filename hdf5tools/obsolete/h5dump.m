function h5dump ( hfile )
% H5DUMP:  metadata dumper
%
% H5DUMP(HFILE) will print metadata for the HDF5 file HFILE.
%
% The amount of indentation for each nested group can be controlled by 
% setting a preference.  To set the group indent to 2 spaces, for 
% example, use
%
%     >> setpref('HDF5TOOLS', 'H5DUMP_INDENT', 2 );
%
% The default value for H5DUMP is 4.
%
% Example:
%     >> h5dump ( 'example.h5' );

%   Copyright 2009 The MathWorks, Inc.




error(nargchk(1,1,nargin,'struct'));

hinfo = hdf5info ( hfile );

options = [];

indent_inc = getpref('MATLAB_IMAGESCI', 'H5DUMP_INDENT', 4);
options.indent_inc = repmat(' ', 1, indent_inc );
options.indent = '';

fprintf ( 1, 'Filename:  %s\n', hinfo.Filename );
dump_group ( hinfo.GroupHierarchy, options );
return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DUMP_GROUP
%     Dumps a single group to the screen.
function dump_group ( group, options ) 

options.indent = [options.indent options.indent_inc];
fprintf ( 1, '%sGroup:  %s\n', options.indent, group.Name );
dump_datatypes ( group.Datatypes, options );
dump_links ( group.Links, options );

options.type = 'group';
dump_attributes ( group.Attributes, options );

options.group_name = group.Name;
dump_datasets ( group.Datasets, options );
dump_groups ( group.Groups, options );
return



function dump_groups ( groups, options )
for j = 1:length(groups)
    dump_group ( groups(j), options );
end
return









%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DUMP_DATASETS
%     Dumps all datasets to the screen.
function dump_datasets ( datasets, options )

if numel(datasets) == 0
    return
end

for j = 1:length(datasets)
    dump_dataset ( datasets(j), options );
end
return








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DUMP_DATASET
%     Dumps a single dataset to the screen.
function dump_dataset ( dataset, options )

options.indent = [options.indent options.indent_inc];

fprintf ( 1, '%sDataset %s\n', options.indent, dataset.Name );

dump_datatype ( dataset.Datatype, options  );
dump_dims ( dataset.Dims, options );
dump_maxdims ( dataset.MaxDims, options );
dump_layout ( dataset.Layout, options );
dump_attributes ( dataset.Attributes, options );
dump_links ( dataset.Links, options );
dump_chunksize ( dataset.Chunksize );

return







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DUMP_LAYOUT
%     Dumps a dataset layout to the screen.
function dump_layout ( layout, options )

options.indent = [options.indent options.indent_inc];

fprintf ( 1, '%sLayout:  %s\n', options.indent, layout );

return







function dump_attributes ( attributes, options )

options.indent = [options.indent options.indent_inc];

if numel(attributes) == 0
    return
end

fprintf ( 1, '%sAttributes:\n', options.indent);
for j = 1:length(attributes)
    dump_attribute(attributes(j), options);
end


return









function shortname = strip_slashes_from_name ( longname )
%
% The name is something like "/Data_Products/VIIRS-Cd-Cov-Agg-IP"
% Strip off the group part.
slashes = findstr(longname,'/');
if isempty(slashes)
	shortname = longname;
else
	shortname = longname(slashes(end)+1:end);
end

return








%===============================================================================
% PARSE_ATTRIBUTE_VALUE
%
% This function creates a string representation of both the attribute type
% and value.  The value will be truncated if it exceeds a certain limit.
function [datatype_string, value_string] = parse_attribute_value ( attval, decimate )


if nargin == 1
    decimate = false;
end

%
% Special case.
if isnumeric ( attval ) && isempty(attval)
    datatype_string = '';
    value_string = '';
    return
end

datatype_string = create_datatype_string(attval);

value_string = create_value_string ( attval, decimate );

if ~ischar(attval) && ( numel(value_string) > 72 )
    [datatype_string, value_string] = parse_attribute_value ( attval, true );
end

return





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% REDUCE_TO_DIMSIZE_ONLY:  
%
% This function reduces an attribute "value" to just an indication of the 
% extent of the dimensions.
function value_string = reduce_to_dimsize_only ( attval )
value_string = sprintf ( '[ %d', size(attval,1) );
for j = 2:ndims(attval)
    value_string = sprintf ( '%s x %d', value_string, size(attval,j) );    
end
value_string = [value_string ' ]' ] ;
return






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DUMP_ATTRIBUTE
%
% Dumps a single attribute to stdout.
function dump_attribute (attribute,options)

options.indent = [options.indent options.indent_inc];

short_name = strip_slashes_from_name ( attribute.Name );
[datatype_string, value_string] = parse_attribute_value ( attribute.Value );

if getpref('MATLAB_IMAGESCI','H5DUMP_VERBOSE',false)
    fprintf ( 1, '%s%s :  %s %s\n', options.indent, short_name, value_string, datatype_string );
else
    fprintf ( 1, '%s%s :  %s\n', options.indent, short_name, value_string );
end
return








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DUMP_DIMS
%     Dumps the dimension extents to the screen.
function dump_dims ( dims, options )

options.indent = [options.indent options.indent_inc];

dimstr = sprintf ( '%d x ', dims );
dimstr = dimstr(1:end-2);
fprintf ( 1, '%sDims:  %s\n', options.indent, dimstr);
return






function dump_maxdims ( dims, options )

options.indent = [options.indent options.indent_inc];

dimstr = sprintf ( '%d x ', dims );
dimstr = dimstr(1:end-2);
fprintf ( 1, '%sMaxDims:  %s\n', options.indent, dimstr);
return






function dump_datatypes ( datatypes, options )

for j = 1:length(datatypes)
    dump_datatype ( datatypes(j), options );
end
return






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DUMP_DATATYPE
%     Dumps the dataset's datatype to the screen.
function dump_datatype ( datatype, options )

options.indent = [options.indent options.indent_inc];

fprintf ( 1, '%sDatatype:  %s\n', options.indent, datatype.Name );
dump_datatype_class ( datatype.Class, options );
dump_datatype_elements ( datatype.Elements, options );

return





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DUMP_DATATYPE_CLASS
%     Dumps the dataset's datatype class to the screen.
function dump_datatype_class ( datatype_class, options )

options.indent = [options.indent options.indent_inc];

if ~isempty(datatype_class)
    fprintf ( 1, '%sClass:  %s\n', options.indent, datatype_class );
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DUMP_DATATYPE_ELEMENT
%     Dumps the dataset's datatype elements to the screen.
function dump_datatype_elements ( datatype_elements, options )

options.indent = [options.indent options.indent_inc];

if ~isempty(datatype_elements)
    fprintf ( 1, '%sFields:  \n%s', options.indent, ...
					elements2str(datatype_elements, options) );
end

return









function dump_chunksize ( chunksize )
if ~isempty(chunksize)
    error ( 'HDF5TOOLS:h5dump:dumpChunkSizeNotImplemented', ...
	        'dump_chunksize not implemented' );
end
return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ELEMENTS2STR
%     Prints out the Elements field of the datatype structure.  
function str = elements2str ( datatype_elements, options )

str = '';
options.indent = [options.indent options.indent_inc];
for j = 1:size(datatype_elements,1)
	str = sprintf('%s%s%s:  %s\n', str, options.indent, datatype_elements{j,1}, datatype_elements{j,2} );
end
return


function dump_links ( links, options )
for j = 1:length(links)
    dump_link ( links(j), options );
end
return



function dump_link ( link, options )
options.indent = [options.indent options.indent_inc];
if link.IsHardLink
    fprintf ( 1, '%s%s is a hardlink to %s\n', options.indent, link.Name, link.Target );
else
    fprintf ( 1, '%s%s is a softlink to %s\n', options.indent, link.Name, link.Target );
end
return



%===============================================================================
% CREATE_DATATYPE_STRING:  Come up with a string representation of the datatype.
function datatype_string = create_datatype_string(attval)

%
% Convert the attribute value into a string representation.
switch class(attval)
case 'hdf5.h5compound'
    datatype_string = 'Compound';
case 'double'
    datatype_string = 'F64';
case 'uint64'
    datatype_string = 'UI64';
case 'int64'
    datatype_string = 'I64';
case 'single'
    datatype_string = 'F32';
case 'int32'
    datatype_string = 'I32';
case 'uint32'
    datatype_string = 'UI32';
case 'int16'
    datatype_string = 'I16';
case 'uint16'
    datatype_string = 'UI16';
case 'int8'
    datatype_string = 'I8';
case 'uint8'
    datatype_string = 'UI8';
case 'hdf5.h5string'
    datatype_string = '';
case 'char'
    datatype_string = 'C';
otherwise
    error ( 'HDF5TOOLS:h5dump:unhandledAttributeClass', ...
	        'attribute value class %s not handled', class(attval ) );
end

return



%===============================================================================
% CREATE_VALUE_STRING:  Creates a string representation of an attribute value.
function value_string = create_value_string ( attval, decimate )
if decimate
       value_string = reduce_to_dimsize_only ( attval );
elseif ~ischar(attval) && (numel(attval) > 1)
       value_string = reduce_to_dimsize_only ( attval );
else
    %
    % Convert the attribute value into a string representation.
    switch class(attval)
	case 'hdf5.h5compound'
    	value_string = '[ compound value ]';
    case 'double'
        value_string = num2str(attval,4);
    case 'single'
        value_string = sprintf ( '%s', num2str(double(attval)) );
    case { 'uint64', 'uint32', 'uint16', 'uint8' }
        value_string = sprintf ( '%s', num2str(double(attval)) );
    case { 'int64', 'int32', 'int16', 'int8' }
        value_string = sprintf ( '%s', num2str(double(attval)) );
    case 'hdf5.h5string'
        value_string = [ '''' attval.Data ''''];
    
    case 'char'
        value_string = attval;
    otherwise
        error ( 'HDF5TOOLS:h5dump:unhandledClass', ...
                'attribute value class ''%s'' not handled', class(attval ) );
    end
end

