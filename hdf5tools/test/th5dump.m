function th5dump()

%   Copyright 2008-2009 The MathWorks, Inc.

fprintf ( 1, '\nStarting H5DUMP tests.\n\n' );
test_example;  
test_dataset_reference;  
test_bitfield;  
fprintf ( 1, '\nFinishing H5DUMP tests.\n\n' );
return


%--------------------------------------------------------------------------
function test_example()

fprintf('\n\n\n');
fprintf('-----------------------------------------------------------------\n');
fprintf('example.h5:  about to dump the example HDF5 file...\n');
pause(2);
fprintf('\n\n\n');

h5dump('example.h5');

fprintf('\n\n\n');
yn = input('example.h5:  Does this look ok (y/n)? ','s');
if strcmp(yn,'n')
	error('Stopping tests at %s', mfilename);
end




%--------------------------------------------------------------------------
function test_dataset_reference()

fprintf('\n\n\n');
fprintf('-----------------------------------------------------------------\n');
fprintf('data/h5ex_t_objref.h5:  about to dump a dataset with references...\n');
pause(2);
fprintf('\n\n\n');

h5dump('data/h5ex_t_objref.h5');

fprintf('\n\n\n');
yn = input('Does this look ok (y/n)? ','s');
if strcmp(yn,'n')
	error('Stopping tests at %s', mfilename);
end




%--------------------------------------------------------------------------
function test_bitfield()

fprintf('\n\n\n');
fprintf('-----------------------------------------------------------------\n');
fprintf('data/h5ex_t_bitf.h5:  about to dump a dataset with bitfields...\n');
pause(2);
fprintf('\n\n\n');

h5dump('data/h5ex_t_bit.h5');

fprintf('\n\n\n');
yn = input('Does this look ok (y/n)? ','s');
if strcmp(yn,'n')
	error('Stopping tests at %s', mfilename);
end
