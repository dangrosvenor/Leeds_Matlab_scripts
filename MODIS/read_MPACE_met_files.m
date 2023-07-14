%run this first (before the read_MPACE_mphys_files script)
savedir='/home/disk/eos1/d.grosvenor/modis_work/plots/';

mpace_dir = '/home/disk/eos1/d.grosvenor/MPACE_ARM_in-situ_cloud_data/arm-iop/poellot-citation/';
%mpace_dir = 'C:\Users\Dan\Documents\logbook\UW\MODIS\MPACE ARM in-situ cloud data\arm-iop\poellot-citation\';
mpace_dir2 = '/home/disk/eos8/d.grosvenor/mpace/';

if exist('box_type')
    switch box_type
        case 'NxN pixel square'
            switch Npix
                case 5
                    load([mpace_dir2 'mpace_flight_mapped_5km_9thOct_2115.mat']);
                    ilinear_mpace = ilinear_mpace_5km;
                case 1

                    %loads the 1km x and y indices with Plat and Plon for the mpace flight
                    %track
                    load([mpace_dir2 'mpace_flight_mapped_1km_9thOct_2115.mat']);
                    ilinear_mpace = ilinear_mpace_1km;
            end
        otherwise
            ilinear_mpace = NaN;
    end
else
    fprintf(1,'\n*** WARNING, no mock L3 data present. ilinear_mpace_xxkm not loaded ***');
end



ncol=49;
NaN_val = 9.99e5;

mpace_met_filename = '04_10_09_20_09_34.mpace';
flight_no='MPACE_09Oct';




eval( [ '[dat_flt' flight_no ',time_flt' flight_no '] = read_mpace_met_files_func(mpace_dir,mpace_met_filename,49,NaN_val);' ] );


mpace_column_numbers

disp('Done read MPACE met files');

