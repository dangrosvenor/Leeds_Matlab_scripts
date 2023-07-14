%N.B. MODIS files are hdf4 not hdf5, so need to use hdfinfo, etc rather
%than hdf5info

%the directory
%filedir='/home/disk/eos10/robwood/';
%filedir='/home/disk/eos1/d.grosvenor/';
%filedir='/home/disk/eos8/d.grosvenor/MAS_data/temp/';
%filedir='/home/disk/eos8/d.grosvenor/';
filedir='/home/disk/eos8/d.grosvenor/MODIS_L3_data/C6/';
filedir='/home/disk/eos8/d.grosvenor/';
filedir='/home/disk/eos5/d.grosvenor/joint_L2/Collection_6/';
filedir='/home/disk/eos15/d.grosvenor/MODIS_L3/aqua_monthly/'; %monthly L3 files from C6.1
filedir='/home/disk/eos15/d.grosvenor/eos8/MOD_L2/Hawaii_Dec2020_ADVANCE/'; %MOD04_L2 files containing AOD data.

%hdf filename
%file_name_h5='MODIS/D3/y2000/MOD08_D3.A2000056.005.2006254074330.hdf';
%file_name_h5='ladsweb.nascom.nasa.gov/allData/5/MOD06_L2/2000/355/MOD06_L2.A2000355.2355.005.2006278230055.hdf';
%file_name_h5='MPACE_JointL2_17days_all_lats/terra/MODATML2.A2004277.0020.051.2010289113442.hdf';
%file_name_h5='ace_mas_980530.017';
%file_name_h5='Arctic_summerL2_2007_20W-60E_70-80N/Joint_5km_files/terra/MODATML2.A2007166.0515.051.2010319070803.hdf';
file_name_h5='y2014_a/MYD08_D3.A2014365.006.2015005181454.hdf';
file_name_h5='MOD_L2/C6/Hamish/MOD06_L2.A2016214.1040.006.2016215022430.hdf'; %Hamish CLARIFY test L2 files
file_name_h5='MOD_L2/C6/Jesus/MYD03.A2014343.1325.006.2014344162340.hdf';
file_name_h5='aqua/2006/001/MYDATML2.A2006001.0450.006.2014224213659.hdf';
file_name_h5='MYD08_M3.A2009060.061.2018042184422.hdf';
file_name_h5='MYD04_L2.A2020356.0055.061.2020357191130.hdf';

%output filename
txt_file='/home/disk/eos1/d.grosvenor/hdf_var_list_L2_Joint_5km.txt';
txt_file='/home/disk/eos1/d.grosvenor/hdf_var_list_L2_AMSRE.txt';
txt_file='/home/disk/eos1/d.grosvenor/idl/LUTs/MODIS_LUT_verylowsun_HDF.txt';
txt_file='/home/disk/eos8/d.grosvenor/primary_prod_Oregon/hdf_file_info.txt';
txt_file='/home/disk/eos8/d.grosvenor/hdf_var_list_L3_Collection_6.txt';
txt_file='/home/disk/eos8/d.grosvenor/hdf_var_list_L2_Collection_6.txt';
txt_file='/home/disk/eos8/d.grosvenor/hdf_var_list_L2_Joint_Collection_6.txt';
%txt_file='/home/disk/eos8/d.grosvenor/hdf_var_list_L2_geolocation.txt';
txt_file=[filedir 'hdf_info_monthly_L3_C6.1.txt'];
txt_file=[filedir 'hdf_info_MYD04_L2_C6.1.txt'];

filename_h5 = [filedir file_name_h5];


%iday=findstr(file_name_h5,'.A');
%modis_year_str=file_name_h5(iday+2:iday+5);
%modis_day_str=file_name_h5(iday+6:iday+8);


SD_id = hdfsd('start',filename_h5,'read'); %open the file
[ndatasets,nglobal_attr,status] = hdfsd('fileinfo',SD_id);

N=ndatasets;


fid=fopen(txt_file,'w');
for nvar=0:N-1
    fprintf(fid,'%d',nvar);
    %    fprintf(fid,' %s',INFO.Vgroup(1).Vgroup(1).SDS(nvar).Name);

    sds_id = hdfsd('select',SD_id, nvar);
    [name,rank,dimsizes,data_type,nattrs,status] = hdfsd('getinfo',sds_id);


    fprintf(fid,' %s',name);
    fprintf(fid,' (dimsizes =');
    for i=1:length(dimsizes)
        fprintf(fid,' %d,',dimsizes(i));
    end
    fprintf(fid,')\n',dimsizes);
    
%    fprintf(fid,' (rank=%d)\n',INFO.Vgroup(1).Vgroup(1).SDS(nvar).Rank);
    %INFO.Vgroup(1).Vgroup(1).SDS(i).Attributes(3).Value; - this contains the
    %detailed name, which might be useful
    %

end
fclose(fid);



disp('Done output hdfinfo');

