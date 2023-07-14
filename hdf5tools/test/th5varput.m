function th5varput()

%   Copyright 2008-2009 The MathWorks, Inc.
test_extend2d();
test_writeFullDouble();
test_writeDoubleWithFillValue();
test_extend1d();
test_writeCompoundVarLengthStrings();

fprintf ( 1, 'All H5VARPUT tests succeeded.\n\n' );
delete('foo.h5');

%--------------------------------------------------------------------------
function test_writeCompoundVarLengthStrings()

v = version('-release');
switch(v)
	case { '2006b', '2007a', '2007b', '2008a', '2008b' }
		fprintf( ['Cannot write compounds with variable length strings ' ...
                  'on 2008b and below.  Filtering ' ...
				  'test_writeCompoundVarLengthStrings test out.\n'] );
		return;

end

copyfile('data/h5ex_t_cmpd.h5','foo.h5');

expData = h5varget('foo.h5','/DS1');
expData.SerialNumber(1) = 999;
h5varput('foo.h5','/DS1',expData);
actData = h5varget('foo.h5','/DS1');

if ~isequal(actData,expData)
	error ( 'failed to verify writing compound dataset');
end

fprintf ( 1, 'Writing a compound dataset succeeded.\n' );

%--------------------------------------------------------------------------
function test_writeDoubleWithFillValue()

create_testfile('foo.h5');

data = rand(8,10);
data(1) = 99;
h5varput('foo.h5','/DS1',data);
outData = h5varget('foo.h5','/DSF');

if ~isnan(outData(1))
	error ( 'failed to verify writing double precision dataset with fill value');
end

fprintf ( 1, 'Writing a double precision variable with fill value succeeded.\n' );

%--------------------------------------------------------------------------
function test_writeFullDouble()
create_testfile('foo.h5');

expVal = rand(8,10);
h5varput('foo.h5','/DS1',expVal);
actVal = h5varget('foo.h5','/DS1');
if ~isequal(actVal,expVal) && isa(actData,'double') 
	error ( 'failed to verify double precision retrieval');
end

fprintf ( 1, 'Writing a double precision variable succeeded.\n' );

%--------------------------------------------------------------------------
function test_extend1d()
create_testfile('foo.h5');

% Ok, the initial size should be 10.  Let's write that, then extend, then
% write some more.
expVal = [0:9];
h5varput('foo.h5','/DSC1',expVal);
actVal = h5varget('foo.h5','/DSC1');
if ~isequal(actVal(:),expVal(:)) && isa(actVal,'double') 
	error ( 'failed to verify double precision 1D retrieval');
end

% Now we write something that will extend the dataset.
expVal = [5:14];
h5varput('foo.h5','/DSC1',5,10,expVal);

fprintf ( 1, 'Writing and extending a double precision 1D succeeded.\n' );

%--------------------------------------------------------------------------
function test_extend2d()
create_testfile('foo.h5');

% Ok, the initial size should be 10.  Let's write that, then extend, then
% write some more.
expVal = ones(10,20);
h5varput('foo.h5','/DSC2',expVal);
actVal = h5varget('foo.h5','/DSC2');
if ~isequal(actVal(:),expVal(:)) && isa(actVal,'double') 
	error ( 'failed to verify double precision 1D retrieval');
end

% Now we write something that will extend the dataset, but only in
% the 2nd dimension.
expVal = 2*ones(10,20);
h5varput('foo.h5','/DSC2',[5 0],[10 20],expVal);
actVal = h5varget('foo.h5','/DSC2',[5 0], [10 20]);

fprintf ( 1, 'Writing and extending a double precision 2D succeeded.\n' );

%--------------------------------------------------------------------------
function create_testfile(hfile)

create_plist = H5P.create('H5P_FILE_CREATE');
file_id = H5F.create(hfile, 'H5F_ACC_TRUNC',create_plist , 'H5P_DEFAULT');

space_id = H5S.create_simple(2,[10 8],[10 8]);

dataset_id = H5D.create(file_id,'/DS1','H5T_NATIVE_DOUBLE', space_id, 'H5P_DEFAULT');
H5D.close(dataset_id);

dtype = H5T.copy('H5T_NATIVE_DOUBLE');
dcpl = H5P.create('H5P_DATASET_CREATE');
H5P.set_fill_value(dcpl,dtype,99);
dataset_id = H5D.create(file_id,'/DSF','H5T_NATIVE_DOUBLE', space_id,dcpl);
H5P.close(dcpl);
H5D.close(dataset_id);
H5T.close(dtype);
H5S.close(space_id);


% Create a 1D dataset with chunking enabled.
dcpl = H5P.create('H5P_DATASET_CREATE');
H5P.set_chunk(dcpl,5);

space_id = H5S.create_simple(1,10,{'H5S_UNLIMITED'});
dtype = H5T.copy('H5T_NATIVE_DOUBLE');
dataset_id = H5D.create(file_id,'/DSC1','H5T_NATIVE_DOUBLE', space_id,dcpl);
H5P.close(dcpl);
H5D.close(dataset_id);
H5T.close(dtype);
H5S.close(space_id);

% Create a 2D dataset with chunking enabled.
dcpl = H5P.create('H5P_DATASET_CREATE');
H5P.set_chunk(dcpl,[5 5]);

space_id = H5S.create_simple(2,[10 20],{'H5S_UNLIMITED','H5S_UNLIMITED'});
dtype = H5T.copy('H5T_NATIVE_DOUBLE');
dataset_id = H5D.create(file_id,'/DSC2','H5T_NATIVE_DOUBLE', space_id,dcpl);
H5P.close(dcpl);
H5D.close(dataset_id);
H5T.close(dtype);
H5S.close(space_id);


H5F.close(file_id);
H5P.close(create_plist);

