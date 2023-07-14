%read_file = '/home/disk/eos1/d.grosvenor/xkcvx_liqCF_004'

%fid=fopen(read_file,'rt'); %is quicker to use fscanf that dlmread
%dat = textscan(fid,'%s %s %s %s %s %s %f %f %f %f %f %f');
%fclose(fid);


%read_file = '/home/disk/eos1/d.grosvenor/xkcvx_liqCF_004'; format_str = '%s %s %s %s %s %s %f %f %f %f %f %f';
  %old output

s='%s ';
f='%f ' ;


clear read_file

%Newer output with Pre-condensation line:-
ifile=0;
ifile=ifile+1; cloud_scheme{ifile}=1; read_file{ifile} = '/home/disk/eos1/d.grosvenor/xkcvg_precond.txt'; format_str{ifile} = [repmat(s,[1 11]) repmat(f,[1 8]) '%s'];
%ifile=ifile+1; cloud_scheme{ifile}=0; read_file{ifile} = '/home/disk/eos1/d.grosvenor/xkcvz_precond.txt'; format_str{ifile} = [repmat(s,[1 11]) repmat(f,[1 8]) '%s'];


%Cloud scheme (1) line:-
ifile=0;
%ifile=ifile+1; cloud_scheme{ifile}=1; read_file{ifile} = '/home/disk/eos1/d.grosvenor/xkcvg_precond.txt'; format_str{ifile} = [repmat(s,[1 11]) repmat(f,[1 8]) '%s'];
ifile=ifile+1; read_file2{ifile} = '/home/disk/eos1/d.grosvenor/xkcvg_cld_scheme01.txt'; format_str2{ifile} = [repmat(s,[1 12]) repmat(f,[1 8])];



%Newer output with Pre-condensation line:-
ifile=0;
%ifile=ifile+1; cloud_scheme{ifile}=0; read_file{ifile} = '/home/disk/eos1/d.grosvenor/UM/output_text/xkcvz002.xkcvz.d14225.t143213.leave_cld_out_precond.txt'; format_str{ifile} = [repmat(s,[1 11]) repmat(f,[1 8]) '%s'];
ifile=ifile+1; cloud_scheme{ifile}=0; read_file{ifile} = '/home/disk/eos1/d.grosvenor/UM/output_text/xkcvh002.xkcvh.d14225.t144039.leave_cld_out_precond.txt'; format_str{ifile} = [repmat(s,[1 11]) repmat(f,[1 8]) '%s'];

%Post-condensation line:-
ifile=0;
%ifile=ifile+1; cloud_scheme{ifile}=0; read_file2{ifile} = '/home/disk/eos1/d.grosvenor/UM/output_text/xkcvz002.xkcvz.d14225.t143213.leave_cld_out_postcond.txt'; format_str2{ifile} = [repmat(s,[1 11]) repmat(f,[1 9]) '%s'];
ifile=ifile+1; cloud_scheme{ifile}=0; read_file2{ifile} = '/home/disk/eos1/d.grosvenor/UM/output_text/xkcvh002.xkcvh.d14225.t144039.leave_cld_out_postcond.txt'; format_str2{ifile} = [repmat(s,[1 11]) repmat(f,[1 9]) '%s'];


%Cloud scheme (1) line:-
ifile=0;
%ifile=ifile+1; read_file3{ifile} = '/home/disk/eos1/d.grosvenor/UM/output_text/xkcvz002.xkcvz.d14225.t143213.leave_cld_out_cld_scheme01.txt'; format_str3{ifile} = [repmat(s,[1 12]) repmat(f,[1 8])];
ifile=ifile+1; read_file3{ifile} = '/home/disk/eos1/d.grosvenor/UM/output_text/xkcvh002.xkcvh.d14225.t144039.leave_cld_out_cld_scheme01.txt'; format_str3{ifile} = [repmat(s,[1 12]) repmat(f,[1 8])];

%Cloud scheme (2) line:-
ifile=0;
%ifile=ifile+1; read_file4{ifile} = '/home/disk/eos1/d.grosvenor/UM/output_text/xkcvz002.xkcvz.d14225.t143213.leave_cld_out_cld_scheme02.txt'; format_str4{ifile} = [repmat(s,[1 10]) repmat(f,[1 6])];
ifile=ifile+1; read_file4{ifile} = '/home/disk/eos1/d.grosvenor/UM/output_text/xkcvh002.xkcvh.d14225.t144039.leave_cld_out_cld_scheme02.txt'; format_str4{ifile} = [repmat(s,[1 10]) repmat(f,[1 6])];





for i=1:length(read_file)
    
    [qV{i},qL{i},qT{i},qs{i},qsL{i},RHtot{i},dmass{i}] = read_cld_scheme_output_func(read_file{i},format_str{i},cloud_scheme{i});   
    [k2{i},qV2{i},qL2{i},qT2{i},qs2{i},dqsdt2{i},qsatfac2{i},dt2{i},qL_new2{i},dqL2{i}] = read_cld_scheme_output_postcond_func(read_file2{i},format_str2{i},cloud_scheme{i});
    [iremove,n{i}, alphaL{i}, aL{i}, bs{i}, Qn{i}, qcl_on_bs{i}, qt{i}, qsl{i}, CF{i}] = read_cld_scheme_output_inside01_func(read_file3{i},format_str3{i}); 
    [qL_init4{i},qL_new4{i},CF4{i},T4{i},TL4{i},T_minus_TL4{i}] = read_cld_scheme_output_inside02_func(read_file4{i},format_str4{i},iremove);
end










%    IF (Qn  <=  -1.0) THEN   !this catch was missing and causing strange values
%       LiqCloudFrac1D(k) = liqCFsmall   !To avoid divide by zero issues.
%       qcl_on_bs= 0.0
%    ELSE IF (Qn  <=  0.) THEN
%       LiqCloudFrac1D(k) = 0.5 * (1. + Qn) * (1. + Qn)
%       qcl_on_bs= (1. + Qn) * (1. + Qn) * (1. + Qn) / 6.
%    ELSE IF (Qn  <   1.) THEN
%       LiqCloudFrac1D(k) = 1. - 0.5 * (1. - Qn) * (1. - Qn)
%       qcl_on_bs=Qn + (1.-Qn) * (1.-Qn) * (1.-Qn)/6.
%    ELSE ! QN  >=  1
%       LiqCloudFrac1D(k) = 1.
%       qcl_on_bs= Qn
%    END IF ! QCN_if